; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_BV, but got a theory atom for that theor
; original: tests/volesti/fuzz_lra/gen/fuzz_0319.smt2 (844 bytes, 14 shrink runs)
; MC volume (original): 3.1956e+07 +- 1.6e+05, ttc: None
;   ttc> c     sec rnd hash      sat     reuse  next         count
;   ttc> Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_BV, but got a theory atom for that theory.
;   ttc> The atom:
;   ttc> (= #b1 __ttc_pbv_0)
; fuzz-meta: {"index": 319, "seed": 1032222, "vars": ["x0", "x1"], "refs": ["x0", "x1"], "box": [["-5000", "3000"], ["-36000", "-11000"]], "features": ["chain_cmp", "comment_tricky", "decl_const", "decl_fun", "define0_bool", "distinct2", "distinct_nary", "exit", "large_consts", "mul_const_left", "neg_unary_const", "num_decimal", "num_int", "num_int_as_decimal", "unary_minus_var"], "mc_volume": 32279999.999999996, "mc_frac": 0.1614}
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-fun x1 () Real)
(define-fun B1 () Bool true)
(assert (and (<= (- 5000) x0) (<= x0 3000) (<= (- 36000) x1) (<= x1 (- 11000))))
