; L-shape [0,2]^2 minus (1,2]x(1,2], written with `not`: area 3.  Non-convex.
; dim: 2
; bbox: -1 3
; exact-volume: 3.0
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (>= x 0) (<= x 2) (>= y 0) (<= y 2)
             (not (and (> x 1) (> y 1)))))
(check-sat)
