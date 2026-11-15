; unit cube, xor of x_i < t_i, t = 0.2,0.4,0.5,0.7,0.9 (2^4 cubes)
; family: parity
; dim: 5
; exact-volume: 0.5
; exact-volume-rational: 1/2
; timeout: 120
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(assert (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1) (>= x3 0) (<= x3 1) (>= x4 0) (<= x4 1)
  (xor (< x0 0.2) (< x1 0.4) (< x2 0.5) (< x3 0.7) (< x4 0.9))))
(check-sat)
