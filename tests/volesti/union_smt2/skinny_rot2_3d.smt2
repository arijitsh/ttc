; rotated crossing slabs 1x1000 extruded by x2 in [0,2]
; family: skinny
; dim: 3
; exact-volume: 3997.9166666666665
; exact-volume-rational: 47975/12
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and (>= x2 0) (<= x2 2)
 (or
  (and (<= (+ (* 3 x0) (* 4 x1)) 2.5) (<= (+ (* (- 3) x0) (* (- 4) x1)) 2.5) (<= (+ (* (- 4) x0) (* 3 x1)) 2500) (<= (+ (* 4 x0) (* (- 3) x1)) 2500))
  (and (<= (+ (* (- 3) x0) (* 4 x1)) 2.5) (<= (+ (* 3 x0) (* (- 4) x1)) 2.5) (<= (+ (* (- 4) x0) (* (- 3) x1)) 2500) (<= (+ (* 4 x0) (* 3 x1)) 2500)))))
(check-sat)
