; class: CRASH
; signature: CRASH: rc=N Symbol 'to_real' not declared as a variable
; original: tests/volesti/fuzz_lra/gen/fuzz_0032.smt2 (882 bytes, 19 shrink runs)
; MC volume (original): 10.0527 +- 0.031, ttc: None
;   ttc> terminate called after throwing an instance of 'cvc5::parser::ParserException'
;   ttc>   what():  Symbol 'to_real' not declared as a variable
; fuzz-meta: {"index": 32, "seed": 1003237, "vars": ["x0", "x1", "x2"], "refs": ["x0", "x1", "x2"], "box": [["3/2", "21/5"], ["-4", "-1/5"], ["4/5", "37/10"]], "features": ["decl_const", "decl_fun", "distinct_nary", "div_var", "eq_nary_neg", "mul_const_right", "multi_assert", "neg_unary_const", "num_decimal", "num_int", "push_declare", "redundant_atom", "to_real", "to_real_in_qflra"], "mc_volume": 10.0776798, "mc_frac": 0.3387}
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-fun x1 () Real)
(declare-fun x2 () Real)
(assert (>= (+ x2 (* x1 (to_real 2))) (- 0.72)))
(declare-const tmpv Real)
(assert (and (<= (/ 3 2) x0) (<= x0 (/ 21 5)) (<= (- 4) x1) (<= x1 (- (/ 1 5))) (<= (/ 4 5) x2) (<= x2 (/ 37 10))))
