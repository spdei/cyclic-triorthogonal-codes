// ====================================================================
// CYCLIC TRIORTHOGONAL CODES IN PRIME DIMENSION p (qudits; also works for p = 2)
// Bron--Kerbosch enumeration of maximal admissible spectral supports,
// constructions A and B, exact distances by Brouwer--Zimmermann.
//
// Adapted from binary-search.m.  Conventions follow the paper:
//   K = spectral support (nonzeros) of the cyclic code, a union of p-cyclotomic cosets mod n.
//   Construction A: stabilizers = cyclic code C with support K, logical row = all-ones.
//       admissible iff  0 notin K, 0 notin 2K, 0 notin 3K
//       i.e.  K meet -K = {}  and  (K+K) meet -K = {}           (same test as the binary script)
//       distance d = min weight of v in Dual(C) with sum(v) ne 0.
//   Construction B: W = cyclic code with support K, punctured at one coordinate.
//       admissible iff  0 notin 3K,  i.e.  (K+K) meet -K = {}    (K meet -K may be nonempty)
//       quantum code [[n-1, 1, d]] with d = d(Dual(W)) - 1.
// No self-orthogonality is imposed (the binary script's self-orthogonality filter is removed).
// ====================================================================

p              := 3;      // prime dimension
n_start        := 4;
n_max          := 121;
n_exact_max    := 121;    // compute exact distances up to this length
max_bk_cosets  := 52;     // above this, fall back to greedy (NOT exhaustive; flagged in the CSV)

F := GF(p);
date_str := Pipe("date '+%Y-%m-%d_%H-%M-%S'", "");
while #date_str gt 0 and date_str[#date_str] eq "\n" do date_str := date_str[1..#date_str-1]; end while;
out_file := Sprintf("Qudit_Triorthogonal_p%o_%o.csv", p, date_str);
fprintf out_file, "p,construction,N,n_quantum,dim,dZ,search,seconds,Roots_K\n";

// ---------------------------------------------------------------- admissibility
function AdmissibleA(n, K)
    if IsEmpty(K) then return true; end if;
    minusK := { (n - x) mod n : x in K };
    if not IsEmpty(K meet minusK) then return false; end if;             // 0 in 2K
    return IsEmpty({ (x + y) mod n : x in K, y in K } meet minusK);      // 0 in 3K
end function;

function AdmissibleB(n, K)
    if IsEmpty(K) then return true; end if;
    minusK := { (n - x) mod n : x in K };
    return IsEmpty({ (x + y) mod n : x in K, y in K } meet minusK);      // 0 in 3K
end function;

function Admissible(cons, n, K)
    if cons eq "A" then return AdmissibleA(n, K); else return AdmissibleB(n, K); end if;
end function;

// ---------------------------------------------------------------- maximal supports
procedure SearchMaximal(cons, R, P, X, cos, n, ~Out)
    if IsEmpty(P) and IsEmpty(X) then
        if not IsEmpty(R) then Append(~Out, R); end if;
        return;
    end if;
    P_iter := P;                       // iterate over a copy: P is modified inside the loop
    for i in P_iter do
        newR := R join cos[i];
        newP := { j : j in P | j ne i and Admissible(cons, n, newR join cos[j]) };
        newX := { j : j in X | Admissible(cons, n, newR join cos[j]) };
        SearchMaximal(cons, newR, newP, newX, cos, n, ~Out);
        Exclude(~P, i); Include(~X, i);
    end for;
end procedure;

function GreedyMaximal(cons, n, cos)
    out := [];
    for seed in cos do
        cur := seed;
        for c in cos do
            if c notsubset cur and Admissible(cons, n, cur join c) then cur join:= c; end if;
        end for;
        if cur notin out then Append(~out, cur); end if;
    end for;
    return out;
end function;

function IsCanonical(n, K)
    s := Sort(Setseq(K));
    for mu in [2..n-1] do
        if Gcd(mu, n) eq 1 and Sort(Setseq({ (x * mu) mod n : x in K })) lt s then return false; end if;
    end for;
    return true;
end function;

// ---------------------------------------------------------------- codes and distances
function CyclicCodeFromSupport(n, K)
    m := Modorder(p, n);
    Fq := GF(p^m);
    Pq<x> := PolynomialRing(Fq);
    // PrimitiveElement, not the field generator: for a prime field (m = 1) the generator is 1.
    beta := PrimitiveElement(Fq)^((p^m - 1) div n);
    g := &*[ x - beta^z : z in {0..n-1} diff K ];
    Pp<y> := PolynomialRing(F);
    return CyclicCode(n, Pp ! [ F ! c : c in Coefficients(g) ]);
end function;

// Construction A: min weight of v in Dual(C) with sum(v) ne 0.
// D = Dual(C + <1>) has codimension 1 in E = Dual(C) (1 notin C because 0 notin K), and
// E \ D is exactly the set of v with sum(v) ne 0.  So d is the least weight w at which E has more
// words than D.  Fast path: the minimum weight of E is attained outside D.
function DistanceA(C, n)
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

// Construction B: d = d(Dual(W)) - 1  (by cyclic symmetry any minimum word can be shifted to the
// punctured coordinate).
function DistanceB(W)
    return MinimumWeight(Dual(W)) - 1;
end function;

// ---------------------------------------------------------------- main loop
for n in [n_start..n_max] do
    if n mod p eq 0 then continue; end if;
    ord := Modorder(p, n);
    unvisited := {1..n-1}; cosets := [];
    while not IsEmpty(unvisited) do
        x0 := Rep(unvisited);
        c := { (x0 * p^i) mod n : i in [0..ord] };
        Append(~cosets, c); unvisited diff:= c;
    end while;

    for cons in ["A", "B"] do
        valid := [ c : c in cosets | Admissible(cons, n, c) ];
        if IsEmpty(valid) then continue; end if;
        maxK := [];
        if #valid le max_bk_cosets then
            mode := "BK";
            SearchMaximal(cons, {}, {1..#valid}, {}, valid, n, ~maxK);
        else
            mode := "GREEDY";
            maxK := GreedyMaximal(cons, n, valid);
        end if;
        for K in maxK do
            if not IsCanonical(n, K) then continue; end if;
            t0 := Cputime();
            C := CyclicCodeFromSupport(n, K);
            assert Dimension(C) eq #K;
            dZ := -1;
            if n le n_exact_max then
                dZ := cons eq "A" select DistanceA(C, n) else DistanceB(C);
            end if;
            nq := cons eq "A" select n else n - 1;
            fprintf out_file, "%o,%o,%o,%o,%o,%o,%o,%o,\"%o\"\n",
                p, cons, n, nq, #K, dZ, mode, Cputime(t0), Sort(Setseq(K));
            printf "p=%o %o n=%o: [[%o,1,%o]] dim %o (%o)\n", p, cons, n, nq, dZ, #K, mode;
        end for;
    end for;
end for;
printf "Done. Results in %o\n", out_file;
