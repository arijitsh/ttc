; A full square plus two zero-volume pieces (a segment x = 3 and a point) that ttc must drop.
; Volume = 1.
; dim: 2
; bbox: -1 4
; exact-volume: 1.0
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (or (and (>= x 0) (<= x 1) (>= y 0) (<= y 1))
            (and (= x 3) (>= y 0) (<= y 2))
            (and (= x 2) (= y 2))))
(check-sat)
