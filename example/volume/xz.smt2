(set-logic QF_LIA)

(declare-const x Int)
(declare-const y Int)

(assert (< x y))
(assert (> x 0))
(assert (< x 5))
(assert (> y 0))
(assert (< y 5))

(check-sat)
(get-model)
