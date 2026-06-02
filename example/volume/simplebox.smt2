(set-logic QF_LIA)
(declare-fun x () Int)
(declare-fun y () Int)
(declare-fun z () Int)

; Derived bounds for x
(assert (>= x -4))
(assert (<= x 5))

; Derived bounds for y
(assert (>= y -2))
(assert (<= y 3))

; Derived bounds for y
(assert (>= z -7))
(assert (<= z 9))

(check-sat)
