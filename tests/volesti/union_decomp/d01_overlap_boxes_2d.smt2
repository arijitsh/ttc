; Two overlapping open squares (0,2)^2 u (1,3)^2: area 4+4-1 = 7.  Overlap of the cover.
; dim: 2
; bbox: -1 4
; exact-volume: 7.0
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (or (and (> x 0) (< x 2) (> y 0) (< y 2))
            (and (> x 1) (< x 3) (> y 1) (< y 3))))
(check-sat)
