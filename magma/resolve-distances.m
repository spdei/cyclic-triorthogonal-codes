// ====================================================================
// RESOLVE DISTANCES THAT THE PYTHON SEARCH LEFT AS LOWER BOUNDS
//
// Codes are given by their generator polynomials over GF(p) (file magma_targets.m), so there is no
// dependence on the choice of primitive root.  Each entry is
//     < p, construction, N, generator coefficients (low to high degree), dim, known distance >
// where, for construction A, "known distance" is exact (SELFTEST) or a rigorous lower bound
// (TARGETS); for construction B the quantum code has length N-1.
//
// Run:   magma resolve-distances.m
// The self-test runs first and aborts if any known distance is not reproduced.
// Results are appended to Resolved_Distances.csv one line at a time, so the run can be stopped
// and resumed (already-resolved entries are skipped on restart).
// ====================================================================

load "magma_targets.m";

out_file := "resolved_distances.csv";
min_p := 3;  max_p := 7;

function CodeFrom(t)
    p := t[1]; N := t[3];
    F := GF(p);
    P<y> := PolynomialRing(F);
    g := P ! [ F ! c : c in t[4] ];
    C := CyclicCode(N, g);
    assert Dimension(C) eq t[5];
    return C, F, N;
end function;

function DistanceA(C, F, n)
    E := Dual(C);
    ones := Matrix(F, 1, n, [ F!1 : i in [1..n] ]);
    D := Dual(LinearCode(VerticalJoin(GeneratorMatrix(C), ones)));
    w := MinimumWeight(E);
    if Dimension(D) eq 0 then return w; end if;
    while w le n do
        if NumberOfWords(E, w) gt NumberOfWords(D, w) then return w; end if;
        w +:= 1;
    end while;
    return -1;
end function;

function Distance(t)
    C, F, N := CodeFrom(t);
    if t[2] eq "A" then return DistanceA(C, F, N); end if;
    return MinimumWeight(Dual(C)) - 1;
end function;

// ---------------------------------------------------------------- self-test
printf "Self-test on %o codes with known exact distance...\n", #SELFTEST;
for t in SELFTEST do
    d := Distance(t);
    nq := t[2] eq "A" select t[3] else t[3] - 1;
    printf "  p=%o %o [[%o,1,%o]]  expected %o  %o\n", t[1], t[2], nq, d, t[6], d eq t[6] select "ok" else "MISMATCH";
    error if d ne t[6], "Self-test failed: the distance routine disagrees with a known value.";
end for;
printf "Self-test passed.\n\n";

// ---------------------------------------------------------------- targets
done := {};
ok := OpenTest(out_file, "r");
if ok then
    for line in Split(Read(out_file), "\n") do
        parts := Split(line, ",");
        if #parts ge 4 and parts[1] ne "p" then Include(~done, parts[1] cat "," cat parts[2] cat "," cat parts[3] cat "," cat parts[#parts]); end if;
    end for;
else
    fprintf out_file, "p,construction,N,n_quantum,dim,previous_lower_bound,dZ_exact,seconds,genkey\n";
end if;

for t in TARGETS do
    if t[1] lt min_p or t[1] gt max_p then continue; end if;
    gk := &cat[ Sprint(c) : c in t[4] ];
    if (Sprintf("%o,%o,%o,", t[1], t[2], t[3]) cat gk) in done then continue; end if;
    nq := t[2] eq "A" select t[3] else t[3] - 1;
    printf "p=%o %o [[%o,1,>=%o]] dim %o ... ", t[1], t[2], nq, t[6], t[5];
    t0 := Cputime();
    d := Distance(t);
    secs := Cputime(t0);
    error if d lt t[6], "Inconsistency: exact distance below the certified lower bound.";
    printf "d = %o  (%o s)\n", d, secs;
    fprintf out_file, "%o,%o,%o,%o,%o,%o,%o,%o,%o\n", t[1], t[2], t[3], nq, t[5], t[6], d, secs, gk;
end for;
printf "Done. Results in %o\n", out_file;
