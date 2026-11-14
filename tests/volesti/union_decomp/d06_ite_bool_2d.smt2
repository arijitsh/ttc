; Boolean ite: if x > 0 then y < x else y < -x, inside [-1,1]^2.  Region = {y < |x|}:
; area = 4 - 1 = 3 (the part y >= |x| is a triangle of area 1).
; dim: 2
; bbox: -2 2
; exact-volume: 3.0
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (<= (- 1) x) (<= x 1) (<= (- 1) y) (<= y 1)
             (ite (> x 0) (< y x) (< y (- x)))))
(check-sat)
