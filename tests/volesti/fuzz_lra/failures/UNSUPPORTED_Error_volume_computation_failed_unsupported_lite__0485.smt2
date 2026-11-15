; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: volume computation failed: unsupported literal encountered during LRA volume conversion
; original: tests/volesti/fuzz_lra/gen/fuzz_0485.smt2 (742 bytes, 12 shrink runs)
; MC volume (original): 2.37534 +- 0.0081, ttc: None
;   ttc> c volume oracle: cg-cdhr
;   ttc> c polytope rounding: yes
;   ttc> c cdd simplification: yes
;   ttc> Error: volume computation failed: unsupported literal encountered during LRA volume conversion
; fuzz-meta: {"index": 485, "seed": 1048988, "vars": ["%m", "x1", "a+b"], "refs": ["%m", "x1", "a+b"], "box": [["-1", "1/2"], ["1", "33/10"], ["14/5", "51/10"]], "features": ["decl_const", "define_args_bool", "mul_const_left", "named", "neg_unary_const", "num_decimal", "num_int", "num_int_as_decimal", "option_after_logic", "special_symbols"], "mc_volume": 2.3935927499999994, "mc_frac": 0.30165}
(declare-const %m Real)
(declare-const x1 Real)
(declare-const a+b Real)
(define-fun inb1 ((u Real) (lo Real) (hi Real)) Bool true)
(assert (and (<= (- 1) %m) (<= %m (/ 1 2)) (<= 1 x1) (<= x1 (/ 33 10)) (<= (/ 14 5) a+b) (<= a+b (/ 51 10))))
