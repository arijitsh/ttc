; class: WRONG
; signature: WRONG: volume zero vs MC (ratio<N)
; original: tests/volesti/fuzz_lra/gen/fuzz_0409.smt2 (741 bytes, 8 shrink runs)
; MC volume (original): 1.10377e+06 +- 3.8e+03, ttc: 0.0
;   ttc> c
;   ttc> c --- [ shutting down ] ----------------------------------------------
;   ttc> c
;   ttc> c exit 0
; fuzz-meta: {"index": 409, "seed": 1041312, "vars": ["x0"], "refs": ["x0"], "box": [["-100000", "3600000"]], "features": ["bbox_not", "bool_ite", "comment_tricky", "const_lhs", "const_moved", "decl_const", "define0_real", "exit", "large_consts", "mul_const_left", "multi_assert", "named", "named_ref", "neg_unary_const", "nested_not", "num_int"], "mc_volume": 1097235.0, "mc_frac": 0.29655}
(set-logic QF_LRA)
(declare-const x0 Real)
(define-fun K1 () Real (- 4240000))
(assert (and (<= (- 100000) x0) (<= x0 3600000)))
