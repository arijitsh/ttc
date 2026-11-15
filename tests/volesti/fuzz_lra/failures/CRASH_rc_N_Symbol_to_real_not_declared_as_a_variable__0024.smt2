; class: CRASH
; signature: CRASH: rc=N Symbol 'to_real' not declared as a variable
; original: tests/volesti/fuzz_lra/gen/fuzz_0024.smt2 (848 bytes, 16 shrink runs)
; MC volume (original): 1.02595 +- 0.0048, ttc: None
;   ttc> terminate called after throwing an instance of 'cvc5::parser::ParserException'
;   ttc>   what():  Symbol 'to_real' not declared as a variable
; fuzz-meta: {"index": 24, "seed": 1002428, "vars": ["x0", "x1", "x2"], "refs": ["x0", "x1", "x2"], "box": [["12/5", "17/5"], ["-16/5", "-2/5"], ["4", "6"]], "features": ["comments", "const_atoms", "decl_const", "eq_nary_neg", "get_model_noopt", "hdr_set_option", "multi_assert", "nary_plus", "neg_unary_const", "num_decimal", "num_int", "num_int_as_decimal", "push_pop", "to_real", "to_real_in_qflra"], "mc_volume": 1.02984, "mc_frac": 0.1839}
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (>= x2 (to_real 4)))
(assert (and (<= (/ 12 5) x0) (<= x0 (/ 17 5)) (<= (- (/ 16 5)) x1) (<= x1 (- (/ 2 5))) (<= 4 x2) (<= x2 6)))
