; class: CRASH
; signature: CRASH: rc=N argument type is not the type of the function's argument type:
; original: tests/volesti/fuzz_lra/gen/fuzz_0369.smt2 (1541 bytes, 42 shrink runs)
; MC volume (original): 0.0028 +- 1.4e-08, ttc: None
;   ttc> argument:  0
;   ttc> has type:  Int
;   ttc> not type: Real
;   ttc> in term : (inb1 x0 (- (/ 7 2500)) 0)
; fuzz-meta: {"index": 369, "seed": 1037272, "vars": ["x0"], "refs": ["x0"], "box": [["-7/2500", "0"]], "features": ["comment_tricky", "comments", "decl_fun", "define_args_bool", "define_args_real", "distinct2", "distinct_nary", "hdr_set_info", "mul_const_left", "mul_const_right", "multi_assert", "neg_unary_const", "no_set_logic", "num_decimal", "num_int", "or", "push_declare", "small_consts"], "mc_volume": 0.0028, "mc_frac": 1.0}
(declare-fun x0 () Real)
(define-fun inb1 ((u Real) (lo Real) (hi Real)) Bool true)
(assert (inb1 x0 (- 0.0028) 0))
(declare-const tmpv Real)
(assert (and (<= (- (/ 7 2500)) x0) (<= x0 0)))
