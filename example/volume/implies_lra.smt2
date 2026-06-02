(set-logic QF_LRA)

(declare-const x Real)
(declare-const y Real)

(assert (and (> x 0) (< x 10)))
(assert (and (> y 0) (< y 10)))
(assert (=> (> x 5) (< y 5)))

(check-sat)
