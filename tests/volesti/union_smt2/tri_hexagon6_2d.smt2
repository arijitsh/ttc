; hexagon (2,0),(1,2),(-1,2),(-2,0),(-1,-2),(1,-2) as 6 triangles from 0
; family: tri
; dim: 2
; exact-volume: 12.0
; exact-volume-rational: 12
; disjuncts: 6
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (<= (* (- 2) x1) 0) (<= (+ (* 2 x0) x1) 4) (<= (+ (* (- 2) x0) x1) 0))
  (and (<= (+ (* 2 x0) (* (- 1) x1)) 0) (<= (* 2 x1) 4) (<= (+ (* (- 2) x0) (* (- 1) x1)) 0))
  (and (<= (+ (* 2 x0) x1) 0) (<= (+ (* (- 2) x0) x1) 4) (<= (* (- 2) x1) 0))
  (and (<= (* 2 x1) 0) (<= (+ (* (- 2) x0) (* (- 1) x1)) 4) (<= (+ (* 2 x0) (* (- 1) x1)) 0))
  (and (<= (+ (* (- 2) x0) x1) 0) (<= (* (- 2) x1) 4) (<= (+ (* 2 x0) x1) 0))
  (and (<= (+ (* (- 2) x0) (* (- 1) x1)) 0) (<= (+ (* 2 x0) (* (- 1) x1)) 4) (<= (* 2 x1) 0))))
(check-sat)
