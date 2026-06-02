(set-logic QF_LRA)
(declare-fun a () Bool)
(declare-fun b () Bool)
(declare-fun x () Real)
(declare-fun y () Real)

(assert (or a b))
(assert (=> a (and (> x 10) (< x 35) (> y 10) (< y 35))))
(assert (=> b (and (> x 15) (< x 40) (> y 15) (< y 40))))

(check-sat)