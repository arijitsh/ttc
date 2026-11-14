; Union of three overlapping triangles (a fan around the origin), rational and decimal
; constants.  T1 = conv{(0,0),(2,0),(0,2)}, T2 = conv{(0,0),(2,0),(2,2)}, T3 = conv{(1,0),(3,0),(2,1)}
; |T1|=|T2|=2, |T1nT2|=1 (tri (0,0),(2,0),(1,1)); |T3|=1, T3n(T1uT2) = T3nT2 = tri (1,0),(2,0),(2,1)
; of area 1/2.  Union area = 3 + 1 - 1/2 = 3.5.
; dim: 2
; bbox: -1 4
; exact-volume: 3.5
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (or (and (>= x 0) (>= y 0) (<= (+ x y) 2))
            (and (<= x 2) (>= y 0) (<= y x))
            (and (>= y 0) (<= y (- x 1)) (<= (+ x y) 3))))
(check-sat)
