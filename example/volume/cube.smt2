(set-logic QF_LIA)
(declare-fun x0 () Int)
(declare-fun x1 () Int)
(declare-fun x2 () Int)

(assert (<= (+ (* 2 x0) (* 3 x1) (* 2 x2)) 10))
(assert (>= (+ (* 2 x0) (* 3 x1)) (- 1)))
(assert (<= (+ (* 4 x0) (* (- 5) x1) (* 2 x2)) 10))
(assert (>= (+ (* 4 x0) (* 5 x2)) 1))

(check-sat)
