; class: UNSUPPORTED
; signature: UNSUPPORTED: Error: volume computation failed: unsupported literal encountered during LRA volume conversion
; original: tests/volesti/fuzz_lra/gen/fuzz_0268.smt2 (730 bytes, 18 shrink runs)
; MC volume (original): 4.9324e+08 +- 7.2e+05, ttc: None
;   ttc> c volume oracle: cg-cdhr
;   ttc> c polytope rounding: yes
;   ttc> c cdd simplification: yes
;   ttc> Error: volume computation failed: unsupported literal encountered during LRA volume conversion
; fuzz-meta: {"index": 268, "seed": 1027071, "vars": ["x0", "x1"], "refs": ["x0", "x1"], "box": [["0", "22000"], ["-5000", "27000"]], "features": ["bbox_chain", "decl_fun", "define_args_bool", "get_value", "hdr_set_option", "large_consts", "mul_const_left", "multi_assert", "neg_unary_const", "num_int", "push_n", "push_unclosed"], "mc_volume": 497481600.0, "mc_frac": 0.70665}
(declare-fun x0 () Real)
(declare-fun x1 () Real)
(define-fun inb1 ((u Real) (lo Real) (hi Real)) Bool true)
(assert (and (<= 0 x0) (<= x0 22000) (<= (- 5000) x1) (<= x1 27000)))
