(set-logic QF_LIA)

(declare-fun b () Bool)
(declare-fun p () Int)

(assert (>= (ite b 10 20) p))
(check-sat)