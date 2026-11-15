; two overlapping triangles (star-of-David-like)
; family: tri
; dim: 2
; exact-volume: 8.0
; exact-volume-rational: 8
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (<= (* (- 4) x1) 0) (<= (+ (* 3 x0) (* 2 x1)) 12) (<= (+ (* (- 3) x0) (* 2 x1)) 0))
  (and (<= (+ (* (- 3) x0) (* (- 2) x1)) (- 4)) (<= (+ (* 3 x0) (* (- 2) x1)) 8) (<= (* 4 x1) 8))))
(check-sat)
