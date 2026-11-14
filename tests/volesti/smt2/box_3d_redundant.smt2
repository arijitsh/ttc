; Unit cube [0,1]^3 plus redundant constraints (duplicates, x+y+z<=10, x<=5, x+y>=-3, x-y<=1).
; dim: 3
; exact-volume: 1.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 0)
  (<= x0 1)
  (>= x1 0)
  (<= x1 1)
  (>= x2 0)
  (<= x2 1)
  (<= x0 1)
  (<= (+ x0 x1 x2) 10)
  (<= x0 5)
  (>= (+ x0 x1) (- 3))
  (<= (- x0 x1) 1)
  (<= (* 2 x2) 2)
))
(check-sat)
