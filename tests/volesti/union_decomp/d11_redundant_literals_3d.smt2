; Box [0,1]^3 described with many redundant literals (weaker bounds, duplicated atoms in
; or/and, a tautological disjunction).  Volume 1.
; dim: 3
; bbox: -1 2
; exact-volume: 1.0
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(declare-const z Real)
(assert (and (>= x 0) (<= x 1) (<= x 2) (<= x 3) (>= x (- 5))
             (>= y 0) (<= y 1) (or (<= y 1) (<= y 0.5)) (or (<= y 2) (>= y 3))
             (>= z 0) (<= z 1) (<= (+ x y z) 3) (<= (+ x y) 2.5)
             (or (<= z 0.5) (> z 0.5))
             (or (< (+ x z) 7) (> x 9))))
(check-sat)
