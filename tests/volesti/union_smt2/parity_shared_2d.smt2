; unit square, (x0<0.3) xor (x0<0.7) xor (x1<0.2): atoms share x0
; family: parity
; dim: 2
; exact-volume: 0.44
; exact-volume-rational: 11/25
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1)
  (xor (< x0 0.3) (< x0 0.7) (< x1 0.2))))
(check-sat)
