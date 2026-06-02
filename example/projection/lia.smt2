(set-logic QF_LIA)

(declare-const p Bool)
(declare-const q Bool)
(declare-const x Int)
(declare-const y Int)
(declare-const k_p Int)
(declare-const k_q Int)

; domain: 0 ≤ x,y < 16  (using only <)
(assert (< -1 x 16))
(assert (< -1 y 16))

; encode (x + y) ≡ 10 (mod 16) and ≡ 5 (mod 16) via fresh ints


(assert (< -1 k_p 2))  ; k_p ∈ {0,1}
(assert (< -1 k_q 2))  ; k_q ∈ {0,1}

(assert (=> p (= (+ x y) (+ 10 (* 16 k_p)))))
(assert (=> q (= (+ x y) (+ 5  (* 16 k_q)))))

(check-sat)
