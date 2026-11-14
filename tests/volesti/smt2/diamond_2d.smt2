; Diamond |x+y|<=1, |x-y|<=1 (square rotated 45 deg), area 2.
; A = [['1', '1'], ['1', '-1']], |det A| = 2; vol = prod(u-l)/|det A|.
; dim: 2
; exact-volume: 2.0
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (and
  (>= (+ x0 x1) (- 1))
  (<= (+ x0 x1) 1)
  (>= (+ x0 (* (- 1) x1)) (- 1))
  (<= (+ x0 (* (- 1) x1)) 1)
))
(check-sat)
