(set-logic QF_LIA)
(declare-fun x () Int)
(declare-fun y () Int)

(assert (and
         (<= (+ (* 2 x) (* 1 y)) 10)
         (>= (+ (* 2 x) (* 1 y)) (- 1))
         (<= (+ (* 4 x) (* (- 5) y)) 10)
         (>= (+ (* 4 x) (* (- 5) y)) 1)
        ))
(check-sat)