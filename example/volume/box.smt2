(set-logic QF_LRA)
(declare-fun x () Real)
(declare-fun y () Real)

(assert (<= (+ (* 2 x) (* 1 y)) 10))
(assert (>= (+ (* 2 x) (* 1 y)) (- 1)))
(assert (<= (+ (* 4 x) (* (- 5) y)) 10))
(assert (>= (+ (* 4 x) (* (- 5) y)) 1))

(check-sat)
