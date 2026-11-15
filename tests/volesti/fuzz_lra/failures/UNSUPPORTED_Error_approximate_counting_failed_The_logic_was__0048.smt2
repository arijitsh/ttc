; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_BV, but got a preprocessing-time term fo
; original: tests/volesti/fuzz_lra/gen/fuzz_0048.smt2 (660 bytes, 11 shrink runs)
; MC volume (original): 0.877811 +- 0.0024, ttc: None
;   ttc> c     sec rnd hash      sat     reuse  next         count
;   ttc> Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_BV, but got a preprocessing-time term for that theory.
;   ttc> The term:
;   ttc> (= #b1 __ttc_pbv_0)
; fuzz-meta: {"index": 48, "seed": 1004851, "vars": ["x0"], "refs": ["x0"], "box": [["-1/2", "17/10"]], "features": ["chain_cmp", "comments", "const_atoms", "decl_fun", "define0_bool", "exit", "mul_const_left", "multi_assert", "neg_unary_const", "num_decimal", "num_int", "push_declare"], "mc_volume": 0.88495, "mc_frac": 0.40225}
(set-logic QF_LRA)
(declare-fun x0 () Real)
(define-fun B1 () Bool (< x0 1.48))
(declare-const tmpv Real)
(assert (and (<= (- (/ 1 2)) x0) (<= x0 (/ 17 10))))
