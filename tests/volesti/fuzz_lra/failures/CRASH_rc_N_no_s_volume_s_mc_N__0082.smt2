; class: CRASH
; signature: CRASH: rc=N no 's volume' (s mc N)
; original: night2/fuzzgen_0/fuzz_0082.smt2 (995 bytes, 36 shrink runs)
; MC volume (original): 1.12051 +- 0.0016, ttc: None
;   ttc> input_stream:7.2: No set-logic command was given before this point.
;   ttc> input_stream:7.2: cvc5 will make all theories available.
;   ttc> input_stream:7.2: Consider setting a stricter logic for (likely) better performance.
;   ttc> input_stream:7.2: To suppress this warning in the future use (set-logic ALL).
; fuzz-meta: {"index": 82, "seed": 1000011282, "vars": ["x0"], "refs": ["x0"], "box": [["8/5", "16/5"]], "features": ["and", "bool_eq", "const_lhs", "decl_fun", "define0_bool", "define_args_real", "hdr_set_info", "let_bool", "mul_const_left", "mul_const_right", "multi_assert", "neg_unary_const", "no_set_logic", "num_decimal", "num_int", "push_pop", "xor"], "mc_volume": 1.11504, "mc_frac": 0.6969}
(declare-fun x0 () Real)
(define-fun B1 () Bool true)
(assert (and (<= (/ 8 5) x0) (<= x0 (/ 16 5))))
