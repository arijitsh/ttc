; L-shape: [0,3]^2 minus (1,3]x(1,3], written with not
; family: shape
; dim: 2
; exact-volume: 5.0
; exact-volume-rational: 5
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and (and (>= x0 0) (<= x0 3) (>= x1 0) (<= x1 3))
  (not (and (> x0 1) (> x1 1)))))
(check-sat)
