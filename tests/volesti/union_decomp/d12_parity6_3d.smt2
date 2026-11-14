; Parity of 6 atoms (two per axis) in the unit cube.  Odd parity per axis on
; x in [.3,.6) (p=.3), y in [.4,.7) (p=.3), z in [.2,.8) (p=.6);
; P(odd total) = (1 - (0.4)(0.4)(-0.2))/2 = 0.516.
; dim: 3
; bbox: -0.5 1.5
; exact-volume: 0.516
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(declare-const z Real)
(assert (and (>= x 0) (<= x 1) (>= y 0) (<= y 1) (>= z 0) (<= z 1)
             (xor (< x 0.3) (< x 0.6) (< y 0.4) (< y 0.7) (< z 0.2) (< z 0.8))))
(check-sat)
