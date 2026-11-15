; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_UF, but got a theory atom for that theor
; original: tests/volesti/fuzz_lra/gen/fuzz_0158.smt2 (968 bytes, 36 shrink runs)
; MC volume (original): 3.7908 +- 0.0069, ttc: None
;   ttc> c     sec rnd hash      sat     reuse  next         count
;   ttc> Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_UF, but got a theory atom for that theory.
;   ttc> The atom:
;   ttc> (= inb1 (lambda ((u Real) (lo Real) (hi Real)) (and (>= (+ u (* (- 1) lo)) 0) (>= (+ (* (- 1) u) hi) 0))))
; fuzz-meta: {"index": 158, "seed": 1015961, "vars": ["x0", "x1"], "refs": ["x0", "x1"], "box": [["1", "9/2"], ["-21/10", "-3/10"]], "features": ["const_lhs", "decl_const", "decl_fun", "define0_bool", "define_args_bool", "exit", "hdr_set_info", "mul_const_left", "mul_const_right", "multi_assert", "nary_plus", "neg_unary_const", "num_decimal", "num_int", "rational_div"], "mc_volume": 3.799215, "mc_frac": 0.60305}
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-fun x1 () Real)
(define-fun inb1 ((u Real) (lo Real) (hi Real)) Bool true)
(define-fun B2 () Bool true)
(assert (and (<= 1 x0) (<= x0 (/ 9 2)) (<= (- (/ 21 10)) x1) (<= x1 (- (/ 3 10)))))
