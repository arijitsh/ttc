(set-logic QF_LRA)

(declare-const x Real)
(declare-const y Real)

(assert (or (and (> x 10) (< x 30) (> y 10) (< y 30))
            (and (> x 50) (< x 70) (> y 50) (< y 70))))

(check-sat)
