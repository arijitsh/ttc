; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_BV, but got a theory atom for that theor
; original: tests/volesti/fuzz_lra/gen/fuzz_0202.smt2 (609 bytes, 15 shrink runs)
; MC volume (original): 0.958872 +- 0.0011, ttc: None
;   ttc> c     sec rnd hash      sat     reuse  next         count
;   ttc> Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_BV, but got a theory atom for that theory.
;   ttc> The atom:
;   ttc> (= #b1 __ttc_pbv_0)
; fuzz-meta: {"index": 202, "seed": 1020405, "vars": ["x0"], "refs": ["x0"], "box": [["-14/5", "-8/5"]], "features": ["comments", "decl_fun", "define0_bool", "mid_check_sat", "mul_const_left", "multi_assert", "neg_unary_const", "nested_not", "no_check_sat", "num_decimal", "num_int_as_decimal", "or"], "mc_volume": 0.9558599999999997, "mc_frac": 0.79655}
(set-logic QF_LRA)
(declare-fun x0 () Real)
(define-fun B1 () Bool true)
(assert (and (<= (- (/ 14 5)) x0) (<= x0 (- (/ 8 5)))))
