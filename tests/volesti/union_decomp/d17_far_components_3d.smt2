; Two disjoint components far apart plus an overlapping third: cube [0,1]^3, cube [50,51]^3
; and the slab [0.5,1.5] x [0,1] x [0,1].  Volume 1 + 1 + 1 - 0.5 = 2.5.
; dim: 3
; bbox: -1 52
; exact-volume: 2.5
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(declare-const z Real)
(assert (or (and (>= x 0) (<= x 1) (>= y 0) (<= y 1) (>= z 0) (<= z 1))
            (and (>= x 50) (<= x 51) (>= y 50) (<= y 51) (>= z 50) (<= z 51))
            (and (>= x 0.5) (<= x 1.5) (>= y 0) (<= y 1) (>= z 0) (<= z 1))))
(check-sat)
