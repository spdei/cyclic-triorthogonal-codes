\\ Cyclotomic cosets and minimal polynomials over F_p, for any prime p and length n with gcd(n,p)=1.
\\ usage: read("cosets_pari_p.gp"); outp(3, 19, "cosets.txt");
outp(p,n,fname)=
{
  my(F, fs, f1, done, cos, s, g, a);
  F = factormod(x^n-1, p)[,1];
  fs = vector(#F, i, F[i]);
  f1 = factormod(polcyclo(n), p)[1,1];
  done = vector(n);
  write(fname, Str("N ", p, " ", n));
  for(r=0, n-1,
    if(!done[r+1],
      cos = List(); s = r;
      while(!done[s+1], done[s+1] = 1; listput(cos, s); s = (p*s) % n);
      if(r == 0,
        g = Mod(1,p)*(x-1),
        a = Mod(Mod(1,p)*x, lift(f1))^r; g = 0;
        for(i=1, #fs, if(subst(lift(fs[i]), x, a) == 0, g = fs[i]; break)));
      write(fname, Str("C ", strjoin(apply(t->Str(t), Vec(cos)), ","), " ",
                       strjoin(apply(t->Str(lift(t)), Vecrev(lift(g))), ",")))));
}
