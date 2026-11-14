; Frame: [-2,2]^2 minus the open square (-1,1)^2: area 16 - 4 = 12.
; dim: 2
; bbox: -3 3
; exact-volume: 12.0
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (<= (- 2) x) (<= x 2) (<= (- 2) y) (<= y 2)
             (not (and (< (- 1) x) (< x 1) (< (- 1) y) (< y 1)))))
(check-sat)
