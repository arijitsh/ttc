; 4 random boxes on a 0.1-grid in [0,1]^6
; family: box
; dim: 6
; exact-volume: 0.028006
; exact-volume-rational: 14003/500000
; disjuncts: 4
; timeout: 120
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(declare-const x5 Real)
(assert (or
  (and (>= x0 0.1) (<= x0 1) (>= x1 0.4) (<= x1 0.6) (>= x2 0) (<= x2 0.4) (>= x3 0.7) (<= x3 1) (>= x4 0.5) (<= x4 0.7) (>= x5 0.4) (<= x5 0.9))
  (and (>= x0 0.3) (<= x0 1) (>= x1 0.6) (<= x1 1) (>= x2 0.1) (<= x2 0.6) (>= x3 0.4) (<= x3 1) (>= x4 0.1) (<= x4 0.9) (>= x5 0.5) (<= x5 0.7))
  (and (>= x0 0.5) (<= x0 1) (>= x1 0.4) (<= x1 0.9) (>= x2 0.1) (<= x2 0.6) (>= x3 0.4) (<= x3 0.6) (>= x4 0) (<= x4 0.5) (>= x5 0.5) (<= x5 1))
  (and (>= x0 0.3) (<= x0 0.9) (>= x1 0) (<= x1 0.7) (>= x2 0.3) (<= x2 0.9) (>= x3 0.6) (<= x3 0.9) (>= x4 0.5) (<= x4 0.7) (>= x5 0.1) (<= x5 0.6))))
(check-sat)
