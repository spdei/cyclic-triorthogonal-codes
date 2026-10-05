// Binary min-weight search with x_0 = 1:  find min |x| s.t. sum_{i in supp x} col_i = target,
// where columns col_i in F_2^K (K <= 1024). Mode A: target = e_{K-1} (last row = all-ones row);
// mode B: target = 0.  Phase 1: complete meet-in-the-middle up to weight wc (rigorous lower bound).
// Phase 2: randomized split MITM for weights > wc (witnesses only). Every hit is verified exactly.
// Output: "L U i1 i2 ..." : distance-search lower bound L (no witness of weight < L), witness weight U
// (or -1) and its support.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define MAXW 16
int n,K,W; uint64_t *col; uint64_t *hc; uint64_t tgt[MAXW*1], htgt;
int mode, oddonly; long captab, caplook; int wmaxr, trials;
typedef struct { uint64_t h, code; } E;
static int cmpE(const void*a,const void*b){ uint64_t x=((E*)a)->h,y=((E*)b)->h; return x<y?-1:x>y; }
E *tab; long ntab, tcap;
int *pool1,np1,*pool2,np2; int h1g,h2g; int I[8],J[8];
uint64_t rnd(){ static uint64_t s=88172645463325252ULL; s^=s<<13; s^=s>>7; s^=s<<17; return s; }
static uint64_t packc(int h,int *A){ uint64_t c=0; for(int j=0;j<h;j++) c|=((uint64_t)A[j])<<(j*10); return c; }
void build(int d,int start,uint64_t h){
  if(d==h1g){ if(ntab==tcap){tcap*=2; tab=realloc(tab,tcap*sizeof(E));} tab[ntab].h=h; tab[ntab].code=packc(d,I); ntab++; return; }
  for(int a=start;a<np1;a++){ I[d]=pool1[a]; build(d+1,a+1,h^hc[pool1[a]]); }
}
int found, wit[64], nwit;
int verify(uint64_t code){
  int A[8]; for(int j=0;j<h1g;j++) A[j]=(code>>(j*10))&1023;
  uint64_t v[MAXW]; memcpy(v,col,sizeof(uint64_t)*W);               // col_0
  for(int j=0;j<h1g;j++) for(int w=0;w<W;w++) v[w]^=col[A[j]*W+w];
  for(int j=0;j<h2g;j++) for(int w=0;w<W;w++) v[w]^=col[J[j]*W+w];
  for(int w=0;w<W;w++) if(v[w]!=tgt[w]) return 0;
  // build support (indices may coincide only if they cancel; handle by xor-count)
  static int cnt[1024]; memset(cnt,0,sizeof(int)*n); cnt[0]^=1;
  for(int j=0;j<h1g;j++) cnt[A[j]]^=1; for(int j=0;j<h2g;j++) cnt[J[j]]^=1;
  nwit=0; for(int i=0;i<n;i++) if(cnt[i]) wit[nwit++]=i;
  return 1;
}
void look(int d,int start,uint64_t h){
  if(found) return;
  if(d==h2g){ E key; key.h=h; E*e=bsearch(&key,tab,ntab,sizeof(E),cmpE); if(!e) return;
    while(e>tab&&(e-1)->h==h) e--; for(;e<tab+ntab&&e->h==h;e++) if(verify(e->code)){found=1;return;} return; }
  for(int a=start;a<np2&&!found;a++){ J[d]=pool2[a]; look(d+1,a+1,h^hc[pool2[a]]); }
}
double C(int a,int b){ if(b<0||b>a) return 0; double r=1; for(int i=0;i<b;i++) r=r*(a-i)/(i+1); return r; }
int attempt(int w,int *P1,int n1,int *P2,int n2){
  int t=w-1; h1g=t/2; h2g=t-h1g; pool1=P1; np1=n1; pool2=P2; np2=n2;
  tcap=1024; ntab=0; tab=malloc(tcap*sizeof(E)); build(0,0,0); qsort(tab,ntab,sizeof(E),cmpE);
  found=0; look(0,0,htgt^hc[0]); free(tab); return found;
}
int main(){
  if(scanf("%d %d %d %ld %ld %d %d",&n,&K,&mode,&captab,&caplook,&wmaxr,&trials)!=7) return 1;
  W=(K+63)/64; col=calloc((size_t)n*W,sizeof(uint64_t)); hc=calloc(n,sizeof(uint64_t));
  uint64_t Rr[1024]; for(int r=0;r<K;r++) Rr[r]=rnd();
  for(int r=0;r<K;r++) for(int i=0;i<n;i++){ int b; if(scanf("%d",&b)!=1) return 1; if(b){ col[i*W+r/64]|=1ULL<<(r%64); hc[i]^=Rr[r]; } }
  memset(tgt,0,sizeof(tgt)); htgt=0;
  if(mode==0){ tgt[(K-1)/64]|=1ULL<<((K-1)%64); htgt=Rr[K-1]; }   // mode A
  oddonly=(mode==0);
  int *all=malloc(sizeof(int)*n); for(int i=1;i<n;i++) all[i-1]=i;
  int L=1;
  for(int w=1; w<=MAXW+20; w++){
    if(oddonly && w%2==0){ continue; }
    int t=w-1, h1=t/2, h2=t-h1;
    if(C(n-1,h1)>captab || C(n-1,h2)>caplook || h1>4 || h2>4) break;
    if(attempt(w,all,n-1,all,n-1)){ printf("%d %d",w,nwit); for(int i=0;i<nwit;i++) printf(" %d",wit[i]); printf("\n"); return 0; }
    L=w+1;
  }
  // Phase 2: randomized split
  int *P1=malloc(sizeof(int)*n),*P2=malloc(sizeof(int)*n);
  for(int w=L; w<=wmaxr; w++){
    if(oddonly && w%2==0) continue;
    int t=w-1, h1=t/2, h2=t-h1; if(h1>4||h2>4) break;
    for(int tr=0;tr<trials;tr++){
      int n1=0,n2=0; for(int i=1;i<n;i++){ if(rnd()&1) P1[n1++]=i; else P2[n2++]=i; }
      if(C(n1,h1)>captab*2 || C(n2,h2)>caplook*2) continue;
      if(attempt(w,P1,n1,P2,n2)){ printf("%d %d",L,nwit); for(int i=0;i<nwit;i++) printf(" %d",wit[i]); printf("\n"); return 0; }
    }
  }
  printf("%d -1\n",L); return 0;
}
