; two parallel rotated slabs 1x1000 overlapping in half their width
; family: skinny
; dim: 2
; exact-volume: 1500.0
; exact-volume-rational: 1500
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(assert (or
  (and (<= (+ (* 5 x0) (* 12 x1)) 6.5) (<= (+ (* (- 5) x0) (* (- 12) x1)) 6.5) (<= (+ (* (- 12) x0) (* 5 x1)) 6500) (<= (+ (* 12 x0) (* (- 5) x1)) 6500))
  (and (<= (+ (* 5 x0) (* 12 x1)) 13) (<= (+ (* (- 5) x0) (* (- 12) x1)) 0) (<= (+ (* (- 12) x0) (* 5 x1)) 6500) (<= (+ (* 12 x0) (* (- 5) x1)) 6500))))
(check-sat)
