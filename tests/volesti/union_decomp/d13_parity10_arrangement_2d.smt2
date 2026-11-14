; Parity of 10 oblique atoms (an arrangement of 10 lines) inside [0,1]^2: many cells / cubes.
; Exact value computed by sampling.
; dim: 2
; bbox: -0.5 1.5
(set-logic QF_LRA)
(declare-const x Real)
(declare-const y Real)
(assert (and (>= x 0) (<= x 1) (>= y 0) (<= y 1)
             (xor (< (+ x y) 1) (< (- x y) 0.2) (< (+ x (* 2 y)) 1.3) (< (- (* 3 x) y) 1)
                  (< (+ (* 0.5 x) y) 0.6) (> (- y x) 0.35) (< x 0.45) (> y 0.15)
                  (< (+ x (* (- 1) (* 0.25 y))) 0.7) (< (+ (* 4 x) (* 3 y)) 2.9))))
(check-sat)
