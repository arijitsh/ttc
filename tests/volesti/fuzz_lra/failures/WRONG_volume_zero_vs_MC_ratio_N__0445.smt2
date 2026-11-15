; class: WRONG
; signature: WRONG: volume zero vs MC (ratio<N)
; original: tests/volesti/fuzz_lra/gen/fuzz_0445.smt2 (766 bytes, 11 shrink runs)
; MC volume (original): 16.6868 +- 0.018, ttc: 0.0
;   ttc> c
;   ttc> c --- [ shutting down ] ----------------------------------------------
;   ttc> c
;   ttc> c exit 0
; fuzz-meta: {"index": 445, "seed": 1044948, "vars": ["x0", "x1", "x2"], "refs": ["x0", "x1", "x2"], "box": [["3/10", "3"], ["-3/2", "6/5"], ["-9/5", "1"]], "features": ["comments", "decl_const", "decl_fun", "define0_real", "exit", "mul_const_left", "nary_plus", "neg_unary_const", "nested_not", "num_decimal", "num_int", "push_pop", "redundant_atom"], "mc_volume": 16.697016, "mc_frac": 0.818}
(declare-const x0 Real)
(declare-const x1 Real)
(declare-fun x2 () Real)
(define-fun K1 () Real 2.51)
(assert (and (<= (/ 3 10) x0) (<= x0 3) (<= (- (/ 3 2)) x1) (<= x1 (/ 6 5)) (<= (- (/ 9 5)) x2) (<= x2 1)))
