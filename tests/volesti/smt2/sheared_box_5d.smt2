; Sheared 5D parallelotope with integer matrix.
; A = [['1', '2', '0', '0', '1'], ['0', '1', '-1', '0', '0'], ['3', '0', '1', '1', '0'], ['0', '0', '2', '1', '-1'], ['1', '1', '1', '1', '2']], |det A| = 18; vol = prod(u-l)/|det A|.
; dim: 5
; exact-volume: 0.8888888888888888
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(declare-const x3 Real)
(declare-const x4 Real)
(assert (and
  (>= (+ x0 (* 2 x1) x4) 0)
  (<= (+ x0 (* 2 x1) x4) 1)
  (>= (+ x1 (* (- 1) x2)) (- 1))
  (<= (+ x1 (* (- 1) x2)) 1)
  (>= (+ (* 3 x0) x2 x3) 2)
  (<= (+ (* 3 x0) x2 x3) 3)
  (>= (+ (* 2 x2) x3 (* (- 1) x4)) 0)
  (<= (+ (* 2 x2) x3 (* (- 1) x4)) 2)
  (>= (+ x0 x1 x2 x3 (* 2 x4)) (- 2))
  (<= (+ x0 x1 x2 x3 (* 2 x4)) 2)
))
(check-sat)
