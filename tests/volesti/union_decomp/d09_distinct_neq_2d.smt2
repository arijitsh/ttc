; Unit square with a measure-zero cut: (distinct x y) and (not (= (+ x y) 1)).  Volume = 1.
; Negated equalities are true almost everywhere; a decomposition must still cover the square.
; dim: 2
; bbox: -1 2
; exact-volume: 1.0
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (>= x 0) (<= x 1) (>= y 0) (<= y 1)
             (distinct x y)
             (not (= (+ x y) 1))))
(check-sat)
