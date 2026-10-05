\\ Binary cyclotomic cosets and minimal polynomials for odd n <= 1023 (PARI/GP >= 2.13 for strjoin).
\\ usage: gp -q cosets_pari.gp
default(parisize,"512M");
out(n)={
  my(F=factormod(x^n-1,2)[,1], fs=[f|f<-F], f1, a, done=vector(n), cos, s, line);
  f1 = factormod(polcyclo(n),2)[1,1];
  write("cosets_p2.txt", Str("N 2 ", n));
  for(r=0, n-1, if(!done[r+1],
     cos=List(); s=r; while(!done[s+1], done[s+1]=1; listput(cos,s); s=(2*s)%n);
     if(r==0, g=Mod(1,2)*(x-1), a=Mod(Mod(1,2)*x, lift(f1))^r; g=0;
        for(i=1,#fs, if(subst(lift(fs[i]),x,a)==0, g=fs[i]; break)));
     my(cs=Vec(lift(g))); cs=Vecrev(cs);
     write("cosets_p2.txt", Str("C ", strjoin(apply(t->Str(t),Vec(cos)),","), " ", strjoin(apply(t->Str(lift(t)),cs),",")))
  ));
}
forstep(n=3,1023,2, out(n));
