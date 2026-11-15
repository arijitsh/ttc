; [0.2,0.5]^3 then [0,1]^3 (small inside big, small listed first)
; family: nest
; dim: 3
; exact-volume: 1.0
; exact-volume-rational: 1
; disjuncts: 2
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (or
  (and (>= x0 0.2) (<= x0 0.5) (>= x1 0.2) (<= x1 0.5) (>= x2 0.2) (<= x2 0.5))
  (and (>= x0 0) (<= x0 1) (>= x1 0) (<= x1 1) (>= x2 0) (<= x2 1))))
(check-sat)
