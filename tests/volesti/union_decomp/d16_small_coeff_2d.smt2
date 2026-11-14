; Atoms with tiny coefficients: 1e-10 x <= 5e-11 means x <= 1/2.  Region [0,1/2] x [0,1]: area 0.5.
; (Probes the absolute kZeroTolerance=1e-9 on the coefficient norm in buildPolytopeMatrices.)
; dim: 2
; bbox: -0.5 1.5
; exact-volume: 0.5
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (>= x 0) (<= x 1) (>= y 0) (<= y 1)
             (<= (* 0.0000000001 x) 0.00000000005)))
(check-sat)
