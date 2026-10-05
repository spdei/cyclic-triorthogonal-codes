// Minimum weight w <= wmax of x in F_p^N with x_0 = 1 and  sum_i x_i b_i = target  (for some listed target),
// where b_i are the columns of B (K x N). Meet-in-the-middle over columns 1..N-1 with EXACT verification:
// every hash hit is reconstructed and compared coordinatewise, so a reported weight comes with a witness.
// Output: "w idx:coef idx:coef ..." or "-1" if none with weight <= wmax.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
int p,K,N,wmax,nt; int *col; int *targets; uint64_t R[512];
typedef struct { uint64_t h, code; } E;
E *tab; long ntab, cap;
static uint64_t hv(const int *v){ uint64_t h=0; for(int i=0;i<K;i++) h+=R[i]*(uint64_t)v[i]; return h; }
static int cmpE(const void*a,const void*b){ uint64_t x=((E*)a)->h,y=((E*)b)->h; return x<y?-1:x>y; }
int acc[40][512]; int idxs[16], cfs[16];
static uint64_t pack(int h){ uint64_t c=0; for(int j=0;j<h;j++){ c|=((uint64_t)idxs[j])<<(j*8); c|=((uint64_t)cfs[j])<<(32+j*4);} return c; }
static void unpack(uint64_t c,int h,int *I,int *C){ for(int j=0;j<h;j++){ I[j]=(c>>(j*8))&255; C[j]=(c>>(32+j*4))&15; } }
int h1g;
void enum1(int d,int start){
  if(d==h1g){ if(ntab==cap){cap*=2; tab=realloc(tab,cap*sizeof(E));} tab[ntab].h=hv(acc[d]); tab[ntab].code=pack(d); ntab++; return; }
  for(int j=start;j<N;j++) for(int c=1;c<p;c++){
    idxs[d]=j; cfs[d]=c; for(int i=0;i<K;i++) acc[d+1][i]=(acc[d][i]+c*col[j*K+i])%p; enum1(d+1,j+1); }
}
int found, h2g, I2[16], C2[16], wI[32], wC[32], wn;
void try_hit(const int *u){
  E key; key.h=hv(u);
  E *e=bsearch(&key,tab,ntab,sizeof(E),cmpE); if(!e) return;
  while(e>tab && (e-1)->h==key.h) e--;
  for(; e<tab+ntab && e->h==key.h; e++){
    int I1[16],C1[16]; unpack(e->code,h1g,I1,C1);
    int v[512]; memset(v,0,sizeof(int)*K);
    for(int j=0;j<h1g;j++) for(int i=0;i<K;i++) v[i]=(v[i]+C1[j]*col[I1[j]*K+i])%p;
    if(memcmp(v,u,sizeof(int)*K)) continue;              // hash collision: reject
    // witness: x_0=1, x_T1 = C1, x_T2 = -C2 (since u = target - b0 + sum_{T2} c col  =>  sum_T1 = u)
    int x[512]; memset(x,0,sizeof(int)*N); x[0]=1;
    for(int j=0;j<h1g;j++) x[I1[j]]=(x[I1[j]]+C1[j])%p;
    for(int j=0;j<h2g;j++) x[I2[j]]=(x[I2[j]]+(p-C2[j]))%p;
    wn=0; for(int i=0;i<N;i++) if(x[i]){ wI[wn]=i; wC[wn]=x[i]; wn++; }
    found=1; return;
  }
}
void enum2(int d,int start){
  if(found) return;
  if(d==h2g){ try_hit(acc[20+d]); return; }
  for(int j=start;j<N && !found;j++) for(int c=1;c<p && !found;c++){
    I2[d]=j; C2[d]=c; for(int i=0;i<K;i++) acc[20+d+1][i]=(acc[20+d][i]+c*col[j*K+i])%p; enum2(d+1,j+1); }
}
int main(){
  if(scanf("%d %d %d %d %d",&p,&K,&N,&wmax,&nt)!=5) return 1;
  col=malloc(sizeof(int)*K*N); targets=malloc(sizeof(int)*K*nt);
  for(int i=0;i<K;i++) for(int j=0;j<N;j++){ int x; if(scanf("%d",&x)!=1) return 1; col[j*K+i]=((x%p)+p)%p; }
  for(int t=0;t<nt;t++) for(int i=0;i<K;i++){ int x; if(scanf("%d",&x)!=1) return 1; targets[t*K+i]=((x%p)+p)%p; }
  srand(987654321); for(int i=0;i<512;i++) R[i]=((uint64_t)rand()<<42)^((uint64_t)rand()<<21)^(uint64_t)rand();
  for(int w=1; w<=wmax; w++){
    int t=w-1; h1g=t/2; h2g=t-h1g;
    cap=1024; ntab=0; tab=malloc(cap*sizeof(E));
    memset(acc[0],0,sizeof(int)*K); enum1(0,1);
    qsort(tab,ntab,sizeof(E),cmpE);
    found=0;
    for(int tt=0; tt<nt && !found; tt++){
      for(int i=0;i<K;i++) acc[20][i]=((targets[tt*K+i]-col[i])%p+p)%p;   // u starts at target - b0
      enum2(0,1);
    }
    free(tab);
    if(found){ printf("%d",wn); for(int i=0;i<wn;i++) printf(" %d:%d",wI[i],wC[i]); printf("\n"); return 0; }
  }
  printf("-1\n"); return 0;
}
