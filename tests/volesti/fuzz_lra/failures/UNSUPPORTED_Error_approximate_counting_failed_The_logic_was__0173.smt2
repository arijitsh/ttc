; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: approximate counting failed: The logic was specified as QF_LIRA, which doesn't include THEORY_BV, but got a theory atom for that theo
; original: tests/volesti/fuzz_lra/gen/fuzz_0173.smt2 (1142 bytes, 9 shrink runs)
; MC volume (original): 8.21272 +- 0.044, ttc: None
;   ttc> c     sec rnd hash      sat     reuse  next         count
;   ttc> Error: approximate counting failed: The logic was specified as QF_LIRA, which doesn't include THEORY_BV, but got a theory atom for that theory.
;   ttc> The atom:
;   ttc> (= #b1 __ttc_pbv_0)
; fuzz-meta: {"index": 173, "seed": 1017476, "vars": ["x0", "x1", "x2", "x3"], "refs": ["x0", "x1", "x2", "x3"], "box": [["-1/2", "9/5"], ["-1/5", "5/2"], ["-39/10", "-1/5"], ["1/2", "29/10"]], "features": ["comments", "const_atoms", "decl_const", "decl_fun", "define0_bool", "div_var", "let_bool", "let_parallel", "mul_const_left", "nary_minus", "nary_plus", "neg_unary_const", "num_decimal", "num_int", "num_int_as_decimal", "to_real", "true_false"], "mc_volume": 8.370980639999997, "mc_frac": 0.1518}
(set-logic QF_LIRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-fun x2 () Real)
(declare-const x3 Real)
(define-fun B1 () Bool true)
(assert (and (<= (- (/ 1 2)) x0) (<= x0 (/ 9 5)) (<= (- (/ 1 5)) x1) (<= x1 (/ 5 2)) (<= (- (/ 39 10)) x2) (<= x2 (- (/ 1 5))) (<= (/ 1 2) x3) (<= x3 (/ 29 10))))
