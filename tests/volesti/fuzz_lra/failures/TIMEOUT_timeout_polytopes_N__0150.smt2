; class: TIMEOUT
; signature: TIMEOUT: timeout (polytopes N)
; original: tests/volesti/fuzz_lra/gen/fuzz_0150.smt2 (1174 bytes, 0 shrink runs)
; MC volume (original): 1.92227 +- 0.0061, ttc: None
;   ttc> c    17.06         8        1.59            18269          19617
;   ttc> c    18.44         9        1.61            17937          19582
;   ttc> c    18.92        10        0.08                0          20476
;   ttc> c    19.46        11        0.09              894          20544
; fuzz-meta: {"index": 150, "seed": 1015154, "vars": ["x0", "x1", "x2"], "refs": ["x0", "x1", "x2"], "box": [["13/5", "19/5"], ["16/5", "11/2"], ["31/10", "26/5"]], "features": ["bbox_chain", "chain_cmp", "decl_fun", "eq_nary_neg", "hdr_set_option", "mul_const_left", "mul_const_right", "multi_assert", "named", "nary_minus", "nary_plus", "neg_unary_const", "no_check_sat", "num_decimal", "num_int", "or", "push_n", "rational_div", "real_ite", "unary_minus_var", "xor"], "mc_volume": 1.9219535999999995, "mc_frac": 0.3316}
(set-option :produce-models true)
(set-logic QF_LRA)
(declare-fun x0 () Real)
(declare-fun x1 () Real)
(declare-fun x2 () Real)
(push 2)
(assert (> x0 8.8))
(pop 2)
(assert (! (xor (< (/ 1437 100) (- (* x1 3) (* x2 (- 2)) x0) 19.29) (or (< (ite (>= 10.97 (- (* x2 (- 1.5)) (* (- 2) x1) (* (- 3) x0)) 8.7) (+ (- x0) (* x1 3) x2 (- 11.57)) (- (* (- 2) x2) (- (/ 218 25)))) 0) (< (ite (<= (* x0 (- 3)) (- (/ 204 25))) (+ (- x2) 4.15) (+ (* x1 3) (- x0) (- 9.4))) 0))) :named n1))
(assert (>= x1 3.2))
(assert (not (= x0 x1 x2 2.6)))
(assert (<= x1 (/ 11 2)))
(assert (and (> (/ 26 5) x2 3.1)
    (>= 3.8 x0 2.6)))
(assert (> (* (- 2) x2) (/ (- 83) 10)))
