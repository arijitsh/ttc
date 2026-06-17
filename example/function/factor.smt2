(set-logic QF_UFBV)

(declare-const w (_ BitVec 12))

; the factorization tool:
; input:  12-bit number
; output: 12-bit pair of factors
(declare-fun mainfunc ((_ BitVec 12)) (_ BitVec 12))

(define-fun out () (_ BitVec 12)
  (mainfunc w))

; extract two 6-bit factors from output
(define-fun x () (_ BitVec 6)
  ((_ extract  5 0) out))

(define-fun y () (_ BitVec 6)
  ((_ extract 11 6) out))

; avoid trivial factors (x and y are 6-bit, so use a 6-bit constant 1)
(assert (bvugt x (_ bv1 6)))
(assert (bvugt y (_ bv1 6)))

; require that the function output is a valid factorization of w
(assert
  (= w
     (bvmul ((_ zero_extend 6) x)
            ((_ zero_extend 6) y))))


(check-sat)
