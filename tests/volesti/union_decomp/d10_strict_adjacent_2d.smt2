; Strict vs non-strict: three boxes that touch along shared faces, one side strict, the other
; non-strict, plus a negated strict literal.  Region = [0,3) x (0,1] u [1,2] x [1,2): area 3 + 1 = 4.
; dim: 2
; bbox: -1 4
; exact-volume: 4.0
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (or (and (<= 0 x) (< x 1) (< 0 y) (<= y 1))
            (and (<= 1 x) (<= x 2) (not (<= y 0)) (<= y 1))
            (and (< 2 x) (< x 3) (< 0 y) (not (> y 1)))
            (and (<= 1 x) (<= x 2) (<= 1 y) (< y 2))))
(check-sat)
