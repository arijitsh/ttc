; let-bound terms, =>, rational constants: s = x/3 + y; (s <= 2/3) => (x - y >= -1/7), in [0,1]^2.
; Region = complement of {s <= 2/3 and x - y < -1/7}.  Exact value computed by sampling.
; dim: 2
; bbox: -0.5 1.5
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (let ((s (+ (* (/ 1 3) x) y)) (inbox (and (>= x 0) (<= x 1) (>= y 0) (<= y 1))))
          (let ((low (<= s (/ 2 3))))
            (and inbox (=> low (>= (- x y) (/ (- 1) 7)))))))
(check-sat)
