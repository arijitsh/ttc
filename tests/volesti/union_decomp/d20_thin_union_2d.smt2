; Two overlapping thin boxes [0,2] x [0,5e-8] u [1,3] x [0,5e-8]: area 3 * 5e-8 = 1.5e-7.
; Anisotropic pieces (aspect 4e7): unit-independence of the decomposition along one axis.
; dim: 2
; bbox: -1 4 -0.00000005 0.0000001
; exact-volume: 1.5e-7
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (>= y 0) (<= y 0.00000005)
             (or (and (>= x 0) (<= x 2)) (and (>= x 1) (<= x 3)))))
(check-sat)
