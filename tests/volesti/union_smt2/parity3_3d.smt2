; unit cube, xor of x_i < t_i, t = 0.5,0.5,0.5 (2^2 cubes)
; family: parity
; dim: 3
; exact-volume: 0.5
; exact-volume-rational: 1/2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1)
  (xor (< x0 0.5) (< x1 0.5) (< x2 0.5))))
(check-sat)
