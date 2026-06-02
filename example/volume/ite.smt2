(set-logic QF_LRA)

(declare-const x Real)
(declare-const y Real)

;; Use ite to express the maximum of x and y equals 10.
(assert (= (ite (> x y) x y) 10))

(check-sat)
