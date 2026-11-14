; Rotated thin rectangle 0<=3x+4y<=50, 0<=-4x+3y<=1 (side 10 x 0.2 rotated, aspect 50).
; A = [['3', '4'], ['-4', '3']], |det A| = 25; vol = prod(u-l)/|det A|.
; dim: 2
; exact-volume: 2.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (>= (+ (* 3 x0) (* 4 x1)) 0)
  (<= (+ (* 3 x0) (* 4 x1)) 50)
  (>= (+ (* (- 4) x0) (* 3 x1)) 0)
  (<= (+ (* (- 4) x0) (* 3 x1)) 1)
))
(check-sat)
