(set-logic QF_LRA)

(declare-const x Real)
(declare-const y Real)

(assert (or (and (> x 0) (< x 10) (> y 0) (< y 10))
            (and (> x 0) (< x 10) (> y 0) (< y 10))
            (and (> x 2) (< x 8) (> y 2) (< y 8))
            (and (> x 50) (< x 60) (> y 50) (< y 60) (<= x 10))))

(check-sat)
