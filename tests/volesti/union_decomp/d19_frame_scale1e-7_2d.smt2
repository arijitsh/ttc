; d03 (frame [-2,2]^2 minus (-1,1)^2) scaled by 1e-7: area 12e-14 = 1.2e-13.
; Decomposition must be unit-independent (Alg. 1 is scale invariant).
; dim: 2
; bbox: -0.0000003 0.0000003
; exact-volume: 1.2e-13
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (<= (- 0.0000002) x) (<= x 0.0000002) (<= (- 0.0000002) y) (<= y 0.0000002)
             (not (and (< (- 0.0000001) x) (< x 0.0000001) (< (- 0.0000001) y) (< y 0.0000001)))))
(check-sat)
