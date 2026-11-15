; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: approximate counting failed: The logic was specified as QF_LIRA, which doesn't include THEORY_BV, but got a preprocessing-time term f
; original: tests/volesti/fuzz_lra/gen/fuzz_0205.smt2 (1593 bytes, 33 shrink runs)
; MC volume (original): 10.4561 +- 0.034, ttc: None
;   ttc> c     sec rnd hash      sat     reuse  next         count
;   ttc> Error: approximate counting failed: The logic was specified as QF_LIRA, which doesn't include THEORY_BV, but got a preprocessing-time term for that theory.
;   ttc> The term:
;   ttc> (= #b1 __ttc_pbv_0)
; fuzz-meta: {"index": 205, "seed": 1020708, "vars": ["x0", "|x)(|", "x2", "|tab\tq|"], "refs": ["x0", "|x)(|", "x2", "|tab\tq|"], "box": [["16/5", "5"], ["23/10", "53/10"], ["-37/10", "-19/10"], ["9/5", "26/5"]], "features": ["bbox_not", "bool_eq", "chain_cmp", "comments", "decl_fun", "define0_bool", "div_var", "dup_assert", "exit", "mul_const_left", "mul_const_right", "multi_assert", "nary_minus", "nary_plus", "neg_unary_const", "nested_not", "num_decimal", "num_int", "or", "push_declare", "push_pop", "quoted_symbols", "rational_div", "to_real"], "mc_volume": 10.322542800000003, "mc_frac": 0.31235}
(set-logic QF_LIRA)
(declare-fun x0 () Real)
(declare-fun |x)(| () Real)
(declare-fun x2 () Real)
(declare-fun |tab	q| () Real)
(define-fun B1 () Bool (<= (+ (* x0 (to_real 3)) (* (- (to_real 2)) |x)(|)) 2.54))
(declare-const tmpv Real)
(assert (and (<= (/ 16 5) x0) (<= x0 5) (<= (/ 23 10) |x)(|) (<= |x)(| (/ 53 10)) (<= (- (/ 37 10)) x2) (<= x2 (- (/ 19 10))) (<= (/ 9 5) |tab	q|) (<= |tab	q| (/ 26 5))))
