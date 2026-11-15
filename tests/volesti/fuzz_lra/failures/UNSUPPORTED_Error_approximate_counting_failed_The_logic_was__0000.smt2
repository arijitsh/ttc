; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: approximate counting failed: The logic was specified as QF_LIRA, which doesn't include THEORY_BV, but got a theory atom for that theo
; original: tests/volesti/fuzz_lra/gen/fuzz_0000.smt2 (1129 bytes, 23 shrink runs)
; MC volume (original): 7.081 +- 0.025, ttc: None
;   ttc> c     sec rnd hash      sat     reuse  next         count
;   ttc> Error: approximate counting failed: The logic was specified as QF_LIRA, which doesn't include THEORY_BV, but got a theory atom for that theory.
;   ttc> The atom:
;   ttc> (= #b1 __ttc_pbv_0)
; fuzz-meta: {"index": 0, "seed": 1000003, "vars": ["x0", "x1", "|x2|", "x3"], "refs": ["|x0|", "|x1|", "x2", "|x3|"], "box": [["-7/5", "12/5"], ["-3/5", "7/10"], ["2/5", "5/2"], ["-7/5", "1"]], "features": ["decl_const", "define0_bool", "define0_real", "distinct2", "exit", "hdr_set_option", "mul_const_left", "nary_plus", "neg_unary_const", "num_decimal", "num_int", "push_unclosed", "quoted_alias", "rational_div", "to_real"], "mc_volume": 7.049755439999999, "mc_frac": 0.28315}
(set-logic QF_LIRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const |x2| Real)
(declare-const x3 Real)
(define-fun B2 () Bool true)
(assert (and (<= (- (/ 7 5)) x0) (<= x0 (/ 12 5)) (<= (- (/ 3 5)) x1) (<= x1 (/ 7 10)) (<= (/ 2 5) |x2|) (<= |x2| (/ 5 2)) (<= (- (/ 7 5)) x3) (<= x3 1)))
