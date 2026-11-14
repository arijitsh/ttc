; Union of the standard simplex {x >= 0, sum <= 1} (vol 1/24), the bounded simplex
; {x_i <= 1/2, sum >= 1} (each x_i >= -1/2; vol 1/24) and the cube [0,1/4]^4 (vol 1/256).
; S1 and S2 meet only on {sum = 1} (measure zero) and the cube lies in S1: union = 1/12.
; dim: 4
; bbox: -1 1.5
; exact-volume: 0.0833333333333
(set-logic QF_LRA)
(declare-const a Real)
(declare-const b Real)
(declare-const c Real)
(declare-const d Real)
(assert (or (and (>= a 0) (>= b 0) (>= c 0) (>= d 0) (<= (+ a b c d) 1))
            (and (<= a 0.5) (<= b 0.5) (<= c 0.5) (<= d 0.5) (>= (+ a b c d) 1))
            (and (>= a 0) (<= a 0.25) (>= b 0) (<= b 0.25) (>= c 0) (<= c 0.25) (>= d 0) (<= d 0.25))))
(check-sat)
