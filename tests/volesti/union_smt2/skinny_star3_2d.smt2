; three crossing slabs 1x300 (aspect 300), directions (1,0),(3,4),(-3,4)
; family: skinny
; dim: 2
; exact-volume: 897.3333333333334
; exact-volume-rational: 2692/3
; disjuncts: 3
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (<= x1 0.5) (<= (* (- 1) x1) 0.5) (<= (* (- 1) x0) 150) (<= x0 150))
  (and (<= (+ (* 4 x0) (* 3 x1)) 2.5) (<= (+ (* (- 4) x0) (* (- 3) x1)) 2.5) (<= (+ (* (- 3) x0) (* 4 x1)) 750) (<= (+ (* 3 x0) (* (- 4) x1)) 750))
  (and (<= (+ (* (- 4) x0) (* 3 x1)) 2.5) (<= (+ (* 4 x0) (* (- 3) x1)) 2.5) (<= (+ (* (- 3) x0) (* (- 4) x1)) 750) (<= (+ (* 3 x0) (* 4 x1)) 750))))
(check-sat)
