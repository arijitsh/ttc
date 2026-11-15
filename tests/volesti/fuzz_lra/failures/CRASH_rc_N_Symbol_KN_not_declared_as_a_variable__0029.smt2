; class: CRASH
; signature: CRASH: rc=N Symbol 'KN' not declared as a variable
; original: tests/volesti/fuzz_lra/gen/fuzz_0029.smt2 (1085 bytes, 28 shrink runs)
; MC volume (original): 1.64536 +- 0.0058, ttc: None
;   ttc> input_stream:2.2: Consider setting a stricter logic for (likely) better performance.
;   ttc> input_stream:2.2: To suppress this warning in the future use (set-logic ALL).
;   ttc> terminate called after throwing an instance of 'cvc5::parser::ParserException'
;   ttc>   what():  Symbol 'K4' not declared as a variable
; fuzz-meta: {"index": 29, "seed": 1002932, "vars": ["|tab\tq|", "x1"], "refs": ["|tab\tq|", "x1"], "box": [["1/2", "14/5"], ["-1", "3/2"]], "features": ["const_atoms", "decl_const", "decl_fun", "define0_real", "dup_assert", "eq_nary", "let", "let_bool", "let_parallel", "mul_const_left", "multi_assert", "neg_unary_const", "no_set_logic", "not", "num_decimal", "num_int", "num_int_as_decimal", "quoted_symbols", "redundant_atom", "true_false"], "mc_volume": 1.65945, "mc_frac": 0.2886}
(declare-const |tab	q| Real)
(declare-fun x1 () Real)
(assert (let ((b2 (> (* (- 1.0) x1) K4))) b1))
(assert (and (<= (/ 1 2) |tab	q|) (<= |tab	q| (/ 14 5)) (<= (- 1) x1) (<= x1 (/ 3 2))))
