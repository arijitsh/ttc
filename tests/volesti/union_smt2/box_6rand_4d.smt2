; 6 random boxes on a 0.1-grid in [0,1]^4
; family: box
; dim: 4
; exact-volume: 0.1504
; exact-volume-rational: 94/625
; disjuncts: 6
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(assert (or
  (and (>= x0 0.3) (<= x0 0.7) (>= x1 0.1) (<= x1 0.9) (>= x2 0.7) (<= x2 0.9) (>= x3 0.1) (<= x3 0.4))
  (and (>= x0 0) (<= x0 0.8) (>= x1 0.4) (<= x1 0.6) (>= x2 0.3) (<= x2 0.9) (>= x3 0.5) (<= x3 0.9))
  (and (>= x0 0.2) (<= x0 1) (>= x1 0.1) (<= x1 0.7) (>= x2 0.3) (<= x2 0.5) (>= x3 0.4) (<= x3 0.8))
  (and (>= x0 0.3) (<= x0 0.6) (>= x1 0.4) (<= x1 0.8) (>= x2 0.5) (<= x2 0.7) (>= x3 0.5) (<= x3 1))
  (and (>= x0 0.3) (<= x0 0.6) (>= x1 0.3) (<= x1 0.8) (>= x2 0.4) (<= x2 0.6) (>= x3 0.4) (<= x3 0.6))
  (and (>= x0 0.4) (<= x0 1) (>= x1 0.4) (<= x1 1) (>= x2 0.3) (<= x2 0.8) (>= x3 0.6) (<= x3 1))))
(check-sat)
