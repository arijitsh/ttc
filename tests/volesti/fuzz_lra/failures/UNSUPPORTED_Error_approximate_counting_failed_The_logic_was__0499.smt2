; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_BV, but got a preprocessing-time term fo
; original: tests/volesti/fuzz_lra/gen/fuzz_0499.smt2 (663 bytes, 10 shrink runs)
; MC volume (original): 0.678657 +- 0.003, ttc: None
;   ttc> c     sec rnd hash      sat     reuse  next         count
;   ttc> Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_BV, but got a preprocessing-time term for that theory.
;   ttc> The term:
;   ttc> (= #b1 __ttc_pbv_0)
; fuzz-meta: {"index": 499, "seed": 1050402, "vars": ["x0"], "refs": ["x0"], "box": [["-19/10", "3/2"]], "features": ["bbox_not", "decl_fun", "define0_bool", "dup_assert", "exit", "hdr_set_option", "mul_const_left", "mul_const_right", "neg_unary_const", "num_decimal", "num_int_as_decimal"], "mc_volume": 0.69054, "mc_frac": 0.2031}
(set-logic QF_LRA)
(declare-fun x0 () Real)
(define-fun B1 () Bool (< (* (- 3.0) x0) (- 1.44)))
(assert (and (<= (- (/ 19 10)) x0) (<= x0 (/ 3 2))))
