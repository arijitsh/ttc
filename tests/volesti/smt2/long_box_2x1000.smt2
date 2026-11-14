; User-reported failure: thin box [0,2] x [0,1000], exact area 2000.
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (>= x 0) (<= x 2) (>= y 0) (<= y 1000)))
(check-sat)
