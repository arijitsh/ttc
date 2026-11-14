; Face x+y+z = 1, x,y,z >= 0 (2D triangle in R^3): 3D volume 0.
; Report-only.
; dim: 3
; exact-volume: 0.0
; report-only: true
(set-logic QF_LRA)
(declare-const x0 Real)
(declare-const x1 Real)
(declare-const x2 Real)
(assert (and
  (>= x0 0)
  (>= x1 0)
  (>= x2 0)
  (= (+ x0 x1 x2) 1)
))
(check-sat)
