; unit square as two triangles split by x0+x1=1
; family: tri
; dim: 2
; exact-volume: 1.0
; exact-volume-rational: 1
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (<= (* (- 1) x1) 0) (<= (+ x0 x1) 1) (<= (* (- 1) x0) 0))
  (and (<= x0 1) (<= x1 1) (<= (+ (* (- 1) x0) (* (- 1) x1)) (- 1)))))
(check-sat)
