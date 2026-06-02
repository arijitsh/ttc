(set-logic QF_LRA)

(declare-const x Real)
(declare-const y Real)

(assert (and (> x 10) (< x 30) (> y 10) (< y 30)))

(check-sat)
