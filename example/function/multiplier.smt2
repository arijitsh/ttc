(set-logic QF_UFBV)

(declare-const x (_ BitVec 6))
(declare-const y (_ BitVec 6))

; unknown implementation of 6-bit multiplication
; input:  x, y are 6-bit
; output: 12-bit product
(declare-fun mainfunc ((_ BitVec 6) (_ BitVec 6)) (_ BitVec 12))

; require the function to match 6x6 -> 12-bit multiplication
(assert
  (= (mainfunc x y)
     (bvmul ((_ zero_extend 6) x)
            ((_ zero_extend 6) y))))

(check-sat)
(get-model)