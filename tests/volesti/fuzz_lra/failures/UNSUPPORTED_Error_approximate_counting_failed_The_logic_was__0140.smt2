; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_UF, but got a theory atom for that theor
; original: tests/volesti/fuzz_lra/gen/fuzz_0140.smt2 (1011 bytes, 22 shrink runs)
; MC volume (original): 1.96125e+12 +- 6.3e+09, ttc: None
;   ttc> c     sec rnd hash      sat     reuse  next         count
;   ttc> Error: approximate counting failed: The logic was specified as QF_LRA, which doesn't include THEORY_UF, but got a theory atom for that theory.
;   ttc> The atom:
;   ttc> (= lf1 (lambda ((|tab	q| Real) (x0 Real)) (+ (* 2 |tab	q|) (* (- 1) x0))))
; fuzz-meta: {"index": 140, "seed": 1014143, "vars": ["x0", "|tab\tq|"], "refs": ["x0", "|tab\tq|"], "box": [["1000000", "3500000"], ["-2600000", "-200000"]], "features": ["comments", "decl_const", "decl_fun", "define0_bool", "define_args_real", "define_param_shadow", "large_consts", "let_parallel", "mul_const_left", "multi_assert", "nary_plus", "neg_unary_const", "no_check_sat", "not", "num_int", "quoted_symbols"], "mc_volume": 1971600000000.0, "mc_frac": 0.3286}
(set-logic QF_LRA)
(declare-fun x0 () Real)
(declare-const |tab	q| Real)
(define-fun lf1 ((|tab	q| Real) (x0 Real)) Real (+ (* (- 1) x0) (* 2 |tab	q|)))
(define-fun B3 () Bool true)
(assert (and (<= 1000000 x0) (<= x0 3500000) (<= (- 2600000) |tab	q|) (<= |tab	q| (- 200000))))
