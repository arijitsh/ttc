; Sheared box 0<=2x+y<=1, 0<=y+z<=2, 0<=x+3z<=3.
; A = [['2', '1', '0'], ['0', '1', '1'], ['1', '0', '3']], |det A| = 7; vol = prod(u-l)/|det A|.
; dim: 3
; exact-volume: 0.8571428571428571
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= (+ (* 2 x0) x1) 0)
  (<= (+ (* 2 x0) x1) 1)
  (>= (+ x1 x2) 0)
  (<= (+ x1 x2) 2)
  (>= (+ x0 (* 3 x2)) 0)
  (<= (+ x0 (* 3 x2)) 3)
))
(check-sat)
