
(set-logic QF_LRA)

(declare-const x Real)
(declare-const y Real)

(assert (or (and (> x 10) (< x 30) (> y 10) (< y 30))
            (and (> x 20) (< x 40) (> y 20) (< y 40))))

(check-sat)
