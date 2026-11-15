; 10 random boxes on a 0.1-grid in [0,1]^3
; family: box
; dim: 3
; exact-volume: 0.669
; exact-volume-rational: 669/1000
; disjuncts: 10
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (or
  (and (>= x0 0.3) (<= x0 0.9) (>= x1 0.2) (<= x1 0.6) (>= x2 0.7) (<= x2 0.9))
  (and (>= x0 0) (<= x0 0.9) (>= x1 0.4) (<= x1 1) (>= x2 0.3) (<= x2 0.6))
  (and (>= x0 0.7) (<= x0 1) (>= x1 0.6) (<= x1 1) (>= x2 0.2) (<= x2 0.5))
  (and (>= x0 0.2) (<= x0 1) (>= x1 0.6) (<= x1 1) (>= x2 0) (<= x2 0.3))
  (and (>= x0 0.2) (<= x0 1) (>= x1 0) (<= x1 0.6) (>= x2 0) (<= x2 0.6))
  (and (>= x0 0.7) (<= x0 1) (>= x1 0.6) (<= x1 0.9) (>= x2 0.7) (<= x2 0.9))
  (and (>= x0 0.5) (<= x0 0.7) (>= x1 0) (<= x1 0.4) (>= x2 0.7) (<= x2 0.9))
  (and (>= x0 0.4) (<= x0 0.9) (>= x1 0.4) (<= x1 0.9) (>= x2 0.6) (<= x2 1))
  (and (>= x0 0.5) (<= x0 1) (>= x1 0.3) (<= x1 0.7) (>= x2 0) (<= x2 0.6))
  (and (>= x0 0.2) (<= x0 0.9) (>= x1 0.5) (<= x1 0.7) (>= x2 0.3) (<= x2 1))))
(check-sat)
