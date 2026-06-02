(set-logic QF_LRA)

(declare-const x Real)
(declare-const y Real)

(assert (and (and (> x 40) (< x 45) (> y 20) (< y 25))
            (and (> x 110) (< x 115) (> y 130) (< y 135))))

(check-sat)
