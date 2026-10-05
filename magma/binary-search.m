// ====================================================================
// BRON-KERBOSCH TRIORTHOGONAL & TRIPLY EVEN CSV DISCOVERY SCRIPT 
// ====================================================================


// Run the shell 'date' command with a clean format
date_str := Pipe("date '+%Y-%m-%d_%H-%M-%S'", "");

// Trim the trailing newline that Pipe includes
while #date_str gt 0 and date_str[#date_str] eq "\n" do
    date_str := date_str[1..#date_str-1];
end while;

// Define your output file
out_file := "Triorthogonal_Codes_Data_" cat date_str cat ".csv";

// Write the CSV Header
fprintf out_file, "n,Dimension,Is_Triply_Even,BCH_Bound,MAX_Bound,Upgrade_Type,c,delta,W,Exact_dZ,A_d,ismax,numNewRows,Roots_K\n";

// Thresholds
n_max := 1023;                // Upper limit for n
n_exact_max := 227;           // Upper limit to compute EXACT quantum distance 
max_bk_cosets := 52; // Max cosets before BK becomes dangerously slow
n_start:= 7;
SetVerbose("Code", 1);



function CalculateQuantumDistance(G0)
    n := NumberOfColumns(G0);
    F2 := GF(2);

    C0     := LinearCode(G0);
    C_perp := Dual(C0);

    // If all-ones is in C0, every dual codeword is orthogonal to it
    // => every dual codeword has even weight => no odd dual word exists
    ones := Vector(F2, [F2!1 : i in [1..n]]);
    if ones in C0 then
        return n + 1;
    end if;

    if Dimension(C_perp) eq 0 then
        return n + 1;
    end if;

    // BZ: global minimum weight of the dual
    d := MinimumWeight(C_perp);

    // Fast path: minimum is already odd
    if d mod 2 eq 1 then
        return d;
    end if;

    // Slow path: minimum is even, search upward for minimum ODD weight.
    // Key insight: Words(C, w : NumWords:=1) uses BZ internally.
    // Since w is odd, ANY word returned has odd weight — no extra parity check needed.
    // We know d_odd >= d+1, and d+1 is odd (since d is even).
    d_odd := d + 1;
    while d_odd le n do
        if #Words(C_perp, d_odd : NumWords := 1) gt 0 then
            return d_odd;
        end if;
        d_odd +:= 2;   // stay on odd values
    end while;

    return n + 1;   // no odd codeword found (should not happen if guard passed)
end function;

function CheckAbsoluteMaximality(G0)
    // Get dimensions and base field
    n := Ncols(G0);
    k := Nrows(G0);
    F := BaseRing(G0); // Assumes GF(2)
    
    // 1. Initialize a sequence for the rows of the constraint matrix M
    // We will build this as a sequence of *sequences* to avoid type errors
    M_rows := [];
    
    // Add the all-ones sequence (enforces even weight for the new row)
    all_ones := [ F!1 : i in [1..n] ];
    Append(~M_rows, all_ones);
    
    // 2. Add all rows of G0 (enforces orthogonality to existing X-stabilizers)
    // Eltseq() converts the matrix row (a vector) into a standard sequence
    for i in [1..k] do
        Append(~M_rows, Eltseq(G0[i]));
    end for;
    
    // 3. Add all pointwise products of pairs of rows in G0 (enforces triorthogonality)
    for i in [1..k-1] do
        for j in [i+1..k] do
            // Build the pointwise product as a standard sequence
            prod := [ G0[i][c] * G0[j][c] : c in [1..n] ];
            Append(~M_rows, prod);
        end for;
    end for;
    
    // 4. Create the Constraint Matrix M (now perfectly type-safe)
    M := Matrix(F, M_rows);
    V := VectorSpace(F, n);
    RowSpace_M := RowSpace(M);
    
    // Valid new rows must be orthogonal to every row in M
    Valid_New_Rows := OrthogonalComplement(V, RowSpace_M);
    Span_G0 := RowSpace(G0);
    
    // 5. Check if the space of valid rows is already entirely contained in G0
    is_maximal := Valid_New_Rows subset Span_G0;
    
    if is_maximal then
        return true, []; // Code cannot be expanded
    else
        // If not maximal, extract the new independent vectors we can add
        new_rows := [];
        
        for v in Basis(Valid_New_Rows) do
            if not (v in Span_G0) then
                Append(~new_rows, v);
                // Dynamically update the span so we only return independent rows
                Span_G0 := sub< V | Span_G0, v >; 
            end if;
        end for;
        
        return false, new_rows;
    end if;
end function;


function CalculateAllBounds(n, K)
    // Roots of the dual code
    R_dual := { (n - x) mod n : x in K };
    
    best_BCH   := 0;
    best_bound := 0;
    best_c     := 1;
    best_L     := 0;   // primary run length (delta in CSV = L+1)
    best_d2    := 0;   // secondary shift count (W in CSV)

    // ----------------------------------------------------------------
    // STEP 1: BCH Bound — longest *consecutive* run in R_dual (step 1)
    // ----------------------------------------------------------------
    for b in R_dual do
        if (b - 1 + n) mod n notin R_dual then   // b starts a run
            L := 0;
            curr := b;
            while curr in R_dual do
                L   +:= 1;
                curr := (curr + 1) mod n;
            end while;
            // L consecutive roots → d ≥ L+1
            if L + 1 gt best_BCH then
                best_BCH := L + 1;
            end if;
        end if;
    end for;

    // ----------------------------------------------------------------
    // STEP 2: Hartmann-Tzeng Bound — rectangle structure
    //
    // For each secondary step s (gcd(s,n)=1):
    //   For each starting root b and primary run length L:
    //     Count d₂ = consecutive shifts j=0,1,...,d₂-1 such that
    //       {b + j*s, b+1 + j*s, ..., b+(L-1) + j*s} ⊆ R_dual
    //     Bound: d ≥ L + d₂
    // ----------------------------------------------------------------
    for s in [1..n-1] do
        if Gcd(s, n) ne 1 then continue; end if;

        for b in R_dual do
            // Only start primary runs at genuine beginnings
            if (b - 1 + n) mod n notin R_dual then

                // Find the primary run length L (consecutive, step 1)
                L := 0;
                curr := b;
                while curr in R_dual do
                    L +:= 1; curr := (curr + 1) mod n;
                end while;

                // Now count how many consecutive secondary shifts
                // shift j: check {b + j*s, ..., b+(L-1) + j*s} ⊆ R_dual
                d2 := 0;
                j  := 0;
                repeat
                    shift_ok := true;
                    for i in [0..L-1] do
                        if (b + i + j*s) mod n notin R_dual then
                            shift_ok := false; break;
                        end if;
                    end for;
                    if shift_ok then d2 +:= 1; j +:= 1;
                    else break;
                    end if;
                until not shift_ok;

                ht_bound := L + d2;   // = (L+1) + (d2-1) in δ/d₂ notation
                if ht_bound gt best_bound then
                    best_bound := ht_bound;
                    best_c     := s;
                    best_L     := L;    // primary run length
                    best_d2    := d2;   // secondary shift count
                end if;
            end if;
        end for;
    end for;

    // ----------------------------------------------------------------
    // STEP 3: Roos Bound — secondary shifts need not be consecutive
    //   For each (b, L, s): instead of requiring shifts 0,1,2,...,d₂-1
    //   to be consecutive, count ALL j ≥ 0 for which the shifted copy
    //   is in R_dual, then take max run of those.
    //   (More complex to implement — flag if found, upgrade = "Roos")
    // ----------------------------------------------------------------
    // [Roos not yet implemented — best_bound already captures true HT]

    max_dist  := Maximum(best_BCH, best_bound);
    upgrade   := "None";
    if max_dist gt best_BCH then
        upgrade := "HT";
    end if;

    // delta = L+1 (= δ in HT notation), W = d₂, c = s
    // Bound formula: delta + W - 1 = (L+1) + d₂ - 1 = L + d₂ ✓
    return max_dist, best_c, best_L + 1, best_d2, upgrade, best_BCH;
end function;

// ---------------------------------------------------------
// FUNCTION: Fast Set-Based Triorthogonality Check
// ---------------------------------------------------------
/* wrong
function IsTriorthogonal(n, K)
    if IsEmpty(K) then return true; end if;
    KplusK := { (x + y) mod n : x in K, y in K };
    minusK := { (-x) mod n : x in K };
    return IsEmpty(KplusK meet minusK);
end function;
*/

function IsValidDefiningSet(n, K)
    if IsEmpty(K) then return true; end if;
    
    // Calculate -K once
    minusK := { (n - x) mod n : x in K };
    
    // 1. Self-Orthogonality Check: K intersect -K must be empty
    if not IsEmpty(K meet minusK) then 
        return false; 
    end if;
    
    // 2. Triorthogonality Check: (K + K) intersect -K must be empty
    KplusK := { (x + y) mod n : x in K, y in K };
    if not IsEmpty(KplusK meet minusK) then 
        return false; 
    end if;
    
    return true;
end function;

// ---------------------------------------------------------
// PROCEDURE: Generalized Bron-Kerbosch
// ---------------------------------------------------------
procedure SearchMaximal(R, P, X, valid_cosets, n, ~MaximalSets)
    if IsEmpty(P) and IsEmpty(X) then
        if not IsEmpty(R) then Append(~MaximalSets, R); end if;
        return;
    end if;
    
    P_iter := P;
    for i in P_iter do
        C := valid_cosets[i];
        new_R := R join C;
        
        new_P := {};
        for j in P do
            if j ne i and IsValidDefiningSet(n, new_R join valid_cosets[j]) then
                Include(~new_P, j);
            end if;
        end for;
        
        new_X := {};
        for j in X do
            if IsValidDefiningSet(n, new_R join valid_cosets[j]) then
                Include(~new_X, j);
            end if;
        end for;
        
        SearchMaximal(new_R, new_P, new_X, valid_cosets, n, ~MaximalSets);
        
        Exclude(~P, i);
        Include(~X, i);
    end for;
end procedure;

// ---------------------------------------------------------
// ENGINE 2: Seed-Based Greedy Search (Fast Fallback)
// ---------------------------------------------------------
function FindGreedyMaximalSets(n, valid_cosets)
    maximal_sets := [];
    for seed in valid_cosets do
        current_set := seed;
        for c in valid_cosets do
            if c ne seed and c notsubset current_set then
                if IsValidDefiningSet(n, current_set join c) then 
                    current_set join:= c; 
                end if;
            end if;
        end for;
        if current_set notin maximal_sets then 
            Append(~maximal_sets, current_set); 
        end if;
    end for;
    return maximal_sets;
end function;

// ---------------------------------------------------------
// FUNCTION: Triply Even Check & Code Construction
// ---------------------------------------------------------
function BuildCodeAndCheckTE(n, K)
    k := #K; 
    if k eq 0 then return false, k, 0; end if; 
    
    m := Modorder(2, n);
    Fq<a_prim> := GF(2^m);
    P<x> := PolynomialRing(Fq);
    prim_root := a_prim^((2^m - 1) div n);
    
    Z := {0..n-1} diff K; 
    g_Fq := &*[ x - prim_root^z : z in Z ];
    
    P2<y> := PolynomialRing(GF(2));
    g_GF2 := P2 ! [ GF(2)!c : c in Coefficients(g_Fq) ];
    
    C := CyclicCode(n, g_GF2);
    G0 := GeneratorMatrix(C);
    B := Basis(C);
    
    is_TE := true;
    for b in B do 
        if Weight(b) mod 8 ne 0 then is_TE := false; break; end if; 
    end for;
    
    if is_TE then
        for i in [1..#B] do
            for j in [i+1..#B] do
                w_int := (Weight(B[i]) + Weight(B[j]) - Weight(B[i] + B[j])) div 2;
                if w_int mod 4 ne 0 then is_TE := false; break; end if;
            end for;
            if not is_TE then break; end if;
        end for;
    end if;

    // Add after the pairwise loop:

    // Replace the entire triple-loop section with this:
    Supps := [ Support(b) : b in B ];
    if is_TE then
        for i in [1..#B] do
            for j in [i+1..#B] do
                for l in [j+1..#B] do
                    // Fast C-level bitset intersection
                    w_triple := #(Supps[i] meet Supps[j] meet Supps[l]);
                    if w_triple mod 2 ne 0 then
                        is_TE := false; break;
                    end if;
                end for;
                if not is_TE then break; end if;
            end for;
            if not is_TE then break; end if;
        end for;
    end if;


    /*if is_TE then
        for i in [1..#B] do
            for j in [i+1..#B] do
                for l in [j+1..#B] do
                    // Triple intersection size
                    w_triple := 0;
                    for c in [1..n] do
                        if B[i][c] eq 1 and B[j][c] eq 1 and B[l][c] eq 1 then
                            w_triple +:= 1;
                        end if;
                    end for;
                    if w_triple mod 2 ne 0 then
                        is_TE := false; break;
                    end if;
                end for;
                if not is_TE then break; end if;
            end for;
            if not is_TE then break; end if;
        end for;
    end if;*/
    
    return is_TE, k, G0;
end function;

// ---------------------------------------------------------
// MAIN EXECUTION
// ---------------------------------------------------------
for n in [n_start..n_max by 2] do
    t0 := Cputime(); 

    ord := Modorder(2, n);
    //if ord gt n/4 and #Factorization(n) lt 2 then continue; end if;

    unvisited := {1..n-1};
    cosets := [];
    while not IsEmpty(unvisited) do
        x := Rep(unvisited);
        C_cos := { (x * (2^i)) mod n : i in [0..ord] };
        Append(~cosets, C_cos);
        unvisited diff:= C_cos;
    end while;
    
    valid_cosets := [ c : c in cosets | IsValidDefiningSet(n, c) ];
    if IsEmpty(valid_cosets) then continue; end if;
    
  printf "Evaluating n = %o | Valid Cosets: %o | ", n, #valid_cosets;
    
    maximal_Ks := [];
    search_type := "";
    
    if #valid_cosets le max_bk_cosets then
        search_type := "BRON-KERBOSCH";
        printf "Mode: %o... ", search_type;
        initial_P := {1..#valid_cosets}; 
        SearchMaximal({}, initial_P, {}, valid_cosets, n, ~maximal_Ks);
    else
        search_type := "GREEDY";
        printf "Mode: %o (Threshold Exceeded)... ", search_type;
        maximal_Ks := FindGreedyMaximalSets(n, valid_cosets);
    end if;


    for K in maximal_Ks do


        // =======================================================
        // DUPLICATE FILTER: Multiplier Equivalence for Sets (Odd n)
        // =======================================================
        is_canonical := true;
        seq_K := Sort(Setseq(K)); // Convert to sorted list for comparison
        
        for mu in [2..n-1] do
            if Gcd(mu, n) eq 1 then
                // Multiply every element in K by mu (modulo n)
                K_mu := { (x * mu) mod n : x in K };
                
                // Compare lexicographically
                seq_K_mu := Sort(Setseq(K_mu));
                if seq_K_mu lt seq_K then
                    is_canonical := false;
                    break;
                end if;
            end if;
        end for;
        
        if not is_canonical then continue; end if; // Skip if it's a duplicate!
        // =======================================================

        // Add after the existing is_canonical check:
        is_self_ortho := IsEmpty(K meet { (n - x) mod n : x in K });
        if not is_self_ortho then continue; end if;

        is_TE, k, G0 := BuildCodeAndCheckTE(n, K);
        if k eq 0 then continue; end if;
        is_max, newrows:=CheckAbsoluteMaximality(G0);
        // Exact Distance
        d_Z := -1; A_d := -1;
        if n le n_exact_max then
            d_Z := CalculateQuantumDistance(G0);
        end if;
        
        // Classical Bounds
        max_dist, step, delta, W, upgrade_str, dist_BCH := CalculateAllBounds(n, K);
        
        K_str := Sprint(K);
        
        fprintf out_file, "%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,\"%o\"\n", 
                n, k, is_TE, dist_BCH, max_dist, upgrade_str, step, delta, W, d_Z, A_d, is_max, #newrows, K_str;
    end for;
    
    t1 := Cputime(t0); 
    printf "Done in %o seconds. Found %o maximal sets.\n", t1, #maximal_Ks;
end for;

printf "\nRun Complete. Data saved to %o.\n", out_file;
