; unit cube, xor of x_i < t_i, t = 0.3,0.6,0.25,0.8 (2^3 cubes)
; family: parity
; dim: 4
; exact-volume: 0.488
; exact-volume-rational: 61/125
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(assert (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1) (>= x3 0) (<= x3 1)
  (xor (< x0 0.3) (< x1 0.6) (< x2 0.25) (< x3 0.8))))
(check-sat)
