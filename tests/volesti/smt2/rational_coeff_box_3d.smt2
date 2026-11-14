; Box written with rational/decimal coefficients: 0<=x/3<=1, 0<=y/4<=1/2, -1<=(5/2)z<=3/2.
; A = [['1/3', '0', '0'], ['0', '1/4', '0'], ['0', '0', '5/2']], |det A| = 5/24; vol = prod(u-l)/|det A|.
; dim: 3
; exact-volume: 6.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= (* (/ 1 3) x0) 0)
  (<= (* (/ 1 3) x0) 1)
  (>= (* (/ 1 4) x1) 0)
  (<= (* (/ 1 4) x1) (/ 1 2))
  (>= (* (/ 5 2) x2) (- 1))
  (<= (* (/ 5 2) x2) (/ 3 2))
))
(check-sat)
