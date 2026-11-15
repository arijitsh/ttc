; class: TIMEOUT
; signature: TIMEOUT: timeout (polytopes N)
; original: tests/volesti/fuzz_lra/gen/fuzz_0328.smt2 (1032 bytes, 0 shrink runs)
; MC volume (original): 18.7401 +- 0.013, ttc: None
;   ttc> c     9.96         2       15.05            26615          26052
;   ttc> c    13.48         3       15.17            26052          26281
;   ttc> c    17.64         4       15.42            26281          26664
;   ttc> c    19.36         5       14.10            18173          16620
; fuzz-meta: {"index": 328, "seed": 1033131, "vars": ["x0", "x1", "x2", "q?r"], "refs": ["x0", "x1", "x2", "q?r"], "box": [["2/5", "11/5"], ["7/2", "59/10"], ["-19/5", "-2/5"], ["39/10", "53/10"]], "features": ["bbox_not", "chain_cmp", "decl_const", "decl_fun", "eq_nary_neg", "hdr_set_option", "implies", "mul_const_left", "named", "nary_plus", "neg_unary_const", "not", "num_decimal", "num_int", "special_symbols", "unary_minus_var"], "mc_volume": 18.658019520000003, "mc_frac": 0.90735}
(set-option :produce-models true)
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-fun x1 () Real)
(declare-const x2 Real)
(declare-const q?r Real)
(assert (and (<= (+ q?r (- x0) (* 2 x2) (- x1)) (- 2.12))
    (>= x1 3.5)
    (<= q?r 5.3)
    (<= x1 5.9)
    (not (< x2 (- 3.8)))
    (not (> x2 (- 0.4)))
    (=> (not (<= (+ (* (- 2) x2) x0) 7.04)) (> (- 21.62) (+ (* (- 3) x1) (* (- 3) x0) (* (- 2) x2) (* (- 3) q?r)) (- 31.06)))
    (>= q?r 3.9)
    (>= x0 0.4)
    (! (not (= x0 x1 x2 q?r 0.4)) :named n1)
    (<= x0 2.2)))
(check-sat)
