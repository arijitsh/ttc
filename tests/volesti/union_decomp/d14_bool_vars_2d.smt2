; Free Boolean variables (projected out): p => x < 1/2, (not p) => y > 1/2, q xor p => x+y < 1.5
; inside [0,1]^2.  sol = exists p,q.  For any x,y take q freely, so region = {x<1/2 or y>1/2}: area 3/4.
; dim: 2
; bbox: -0.5 1.5
; exact-volume: 0.75
; report-only: true
; (ttc routes formulas with free Bool variables to Boolean model counting, not to volume.)
(set-logic QF_LRA)
(declare-const p Bool)
(declare-const q Bool)
(declare-const x Real)
(declare-const y Real)
(assert (and (>= x 0) (<= x 1) (>= y 0) (<= y 1)
             (=> p (< x 0.5)) (=> (not p) (> y 0.5))
             (=> (xor q p) (< (+ x y) 1.5))))
(check-sat)
