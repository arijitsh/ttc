(set-logic QF_LRA)

; Booleans
(declare-fun a () Bool)
(declare-fun b () Bool)
(declare-fun c () Bool)

; Reals
(declare-fun A () Real)
(declare-fun B () Real)
(declare-fun C () Real)
(declare-fun X () Real)
(declare-fun Y () Real)
(declare-fun Z () Real)

; X = a*A + c*C
(assert (= X (+ (ite a A 0.0) (ite c C 0.0))))

; Y = -b*B - a*A
(assert (= Y (+ (ite b (- B) 0.0) (ite a (- A) 0.0))))

; Z = b*B - c*C
(assert (= Z (+ (ite b B 0.0) (ite c (- C) 0.0))))

(assert (<= X 1.0))
(assert (<= Y 1.0))
(assert (<= Z 1.0))

(check-sat)
(get-model)