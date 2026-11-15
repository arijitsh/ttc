; two crossing slabs 1x100 along (4,-3)/5 and (4,3)/5
; family: skinny
; dim: 2
; exact-volume: 198.95833333333334
; exact-volume-rational: 4775/24
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (<= (+ (* 3 x0) (* 4 x1)) 2.5) (<= (+ (* (- 3) x0) (* (- 4) x1)) 2.5) (<= (+ (* (- 4) x0) (* 3 x1)) 250) (<= (+ (* 4 x0) (* (- 3) x1)) 250))
  (and (<= (+ (* (- 3) x0) (* 4 x1)) 2.5) (<= (+ (* 3 x0) (* (- 4) x1)) 2.5) (<= (+ (* (- 4) x0) (* (- 3) x1)) 250) (<= (+ (* 4 x0) (* 3 x1)) 250))))
(check-sat)
