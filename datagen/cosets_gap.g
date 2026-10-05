# Cyclotomic cosets and minimal polynomials over GF(p) for lengths n coprime to p.
# usage (in GAP):  Read("cosets_gap.g"); WriteCosets(5, [4,6,7,8,9,11], "cosets_p5.txt");
DoN := function(p, N, out)
  local m, q, a, cos, done, i, c, mp;
  m := OrderMod(p, N); q := p^m; a := Z(q)^((q-1)/N);
  done := []; AppendTo(out, "N ", p, " ", N, "\n");
  for i in [0..N-1] do
    if not i in done then
      c := Set(List([0..m-1], j -> (i*p^j) mod N)); Append(done, c);
      mp := MinimalPolynomial(GF(p), a^i);
      AppendTo(out, "C ", JoinStringsWithSeparator(List(c, String), ","), " ",
        JoinStringsWithSeparator(List(CoefficientsOfUnivariatePolynomial(mp), x -> String(IntFFE(x))), ","), "\n");
    fi;
  od;
end;
WriteCosets := function(p, lengths, fname)
  local out, N;
  out := OutputTextFile(fname, false); SetPrintFormattingStatus(out, false);
  for N in lengths do DoN(p, N, out); od;
  CloseStream(out);
end;
