;; SMT-LRA encoding for power-flow toy model: node14_edge20_gen100_numgen30_flowlim50_seed0
;; Constants inlined; one flow variable per undirected line; f_j_i is -f_i_j.
;; Nodes: 14; Edges: 20

(set-logic QF_LRA)

;; ---------- Edge variables (canonical orientation i>j) ----------
;; Edge {11,4} with capacity 48
(declare-fun proj_x_11_4 () Bool)
(declare-fun f_11_4 () Real)
;; Edge {12,11} with capacity 50
(declare-fun proj_x_12_11 () Bool)
(declare-fun f_12_11 () Real)
;; Edge {11,7} with capacity 47
(declare-fun proj_x_11_7 () Bool)
(declare-fun f_11_7 () Real)
;; Edge {11,6} with capacity 47
(declare-fun proj_x_11_6 () Bool)
(declare-fun f_11_6 () Real)
;; Edge {6,1} with capacity 51
(declare-fun proj_x_6_1 () Bool)
(declare-fun f_6_1 () Real)
;; Edge {8,7} with capacity 49
(declare-fun proj_x_8_7 () Bool)
(declare-fun f_8_7 () Real)
;; Edge {14,7} with capacity 52
(declare-fun proj_x_14_7 () Bool)
(declare-fun f_14_7 () Real)
;; Edge {9,6} with capacity 51
(declare-fun proj_x_9_6 () Bool)
(declare-fun f_9_6 () Real)
;; Edge {5,4} with capacity 46
(declare-fun proj_x_5_4 () Bool)
(declare-fun f_5_4 () Real)
;; Edge {10,9} with capacity 49
(declare-fun proj_x_10_9 () Bool)
(declare-fun f_10_9 () Real)
;; Edge {4,2} with capacity 48
(declare-fun proj_x_4_2 () Bool)
(declare-fun f_4_2 () Real)
;; Edge {13,11} with capacity 49
(declare-fun proj_x_13_11 () Bool)
(declare-fun f_13_11 () Real)
;; Edge {3,2} with capacity 50
(declare-fun proj_x_3_2 () Bool)
(declare-fun f_3_2 () Real)
;; Edge {8,2} with capacity 52
(declare-fun proj_x_8_2 () Bool)
(declare-fun f_8_2 () Real)
;; Edge {8,3} with capacity 47
(declare-fun proj_x_8_3 () Bool)
(declare-fun f_8_3 () Real)
;; Edge {9,2} with capacity 47
(declare-fun proj_x_9_2 () Bool)
(declare-fun f_9_2 () Real)
;; Edge {11,8} with capacity 49
(declare-fun proj_x_11_8 () Bool)
(declare-fun f_11_8 () Real)
;; Edge {5,2} with capacity 48
(declare-fun proj_x_5_2 () Bool)
(declare-fun f_5_2 () Real)
;; Edge {11,2} with capacity 48
(declare-fun proj_x_11_2 () Bool)
(declare-fun f_11_2 () Real)
;; Edge {3,1} with capacity 51
(declare-fun proj_x_3_1 () Bool)
(declare-fun f_3_1 () Real)

;; ---------- Line-outage and capacity limits ----------
(assert (=> (not proj_x_11_4) (= f_11_4 0.0)))
(assert (and (<= (- 48.0) f_11_4) (<= f_11_4 48.0)))
(assert (=> (not proj_x_12_11) (= f_12_11 0.0)))
(assert (and (<= (- 50.0) f_12_11) (<= f_12_11 50.0)))
(assert (=> (not proj_x_11_7) (= f_11_7 0.0)))
(assert (and (<= (- 47.0) f_11_7) (<= f_11_7 47.0)))
(assert (=> (not proj_x_11_6) (= f_11_6 0.0)))
(assert (and (<= (- 47.0) f_11_6) (<= f_11_6 47.0)))
(assert (=> (not proj_x_6_1) (= f_6_1 0.0)))
(assert (and (<= (- 51.0) f_6_1) (<= f_6_1 51.0)))
(assert (=> (not proj_x_8_7) (= f_8_7 0.0)))
(assert (and (<= (- 49.0) f_8_7) (<= f_8_7 49.0)))
(assert (=> (not proj_x_14_7) (= f_14_7 0.0)))
(assert (and (<= (- 52.0) f_14_7) (<= f_14_7 52.0)))
(assert (=> (not proj_x_9_6) (= f_9_6 0.0)))
(assert (and (<= (- 51.0) f_9_6) (<= f_9_6 51.0)))
(assert (=> (not proj_x_5_4) (= f_5_4 0.0)))
(assert (and (<= (- 46.0) f_5_4) (<= f_5_4 46.0)))
(assert (=> (not proj_x_10_9) (= f_10_9 0.0)))
(assert (and (<= (- 49.0) f_10_9) (<= f_10_9 49.0)))
(assert (=> (not proj_x_4_2) (= f_4_2 0.0)))
(assert (and (<= (- 48.0) f_4_2) (<= f_4_2 48.0)))
(assert (=> (not proj_x_13_11) (= f_13_11 0.0)))
(assert (and (<= (- 49.0) f_13_11) (<= f_13_11 49.0)))
(assert (=> (not proj_x_3_2) (= f_3_2 0.0)))
(assert (and (<= (- 50.0) f_3_2) (<= f_3_2 50.0)))
(assert (=> (not proj_x_8_2) (= f_8_2 0.0)))
(assert (and (<= (- 52.0) f_8_2) (<= f_8_2 52.0)))
(assert (=> (not proj_x_8_3) (= f_8_3 0.0)))
(assert (and (<= (- 47.0) f_8_3) (<= f_8_3 47.0)))
(assert (=> (not proj_x_9_2) (= f_9_2 0.0)))
(assert (and (<= (- 47.0) f_9_2) (<= f_9_2 47.0)))
(assert (=> (not proj_x_11_8) (= f_11_8 0.0)))
(assert (and (<= (- 49.0) f_11_8) (<= f_11_8 49.0)))
(assert (=> (not proj_x_5_2) (= f_5_2 0.0)))
(assert (and (<= (- 48.0) f_5_2) (<= f_5_2 48.0)))
(assert (=> (not proj_x_11_2) (= f_11_2 0.0)))
(assert (and (<= (- 48.0) f_11_2) (<= f_11_2 48.0)))
(assert (=> (not proj_x_3_1) (= f_3_1 0.0)))
(assert (and (<= (- 51.0) f_3_1) (<= f_3_1 51.0)))

;; ---------- Nodal balance ----------
(assert (= (+ (- f_6_1) (- f_3_1)) 34.0))  ;; node 1
(assert (= (+ (- f_4_2) (- f_3_2) (- f_8_2) (- f_9_2) (- f_5_2) (- f_11_2)) (- 18.0)))  ;; node 2
(assert (= (+ f_3_2 (- f_8_3) f_3_1) (- 10.0)))  ;; node 3
(assert (= (+ (- f_11_4) (- f_5_4) f_4_2) (- 9.0)))  ;; node 4
(assert (= (+ f_5_4 f_5_2) (- 2.0)))  ;; node 5
(assert (= (+ (- f_11_6) f_6_1 (- f_9_6)) (- 7.0)))  ;; node 6
(assert (= (+ (- f_11_7) (- f_8_7) (- f_14_7)) 29.0))  ;; node 7
(assert (= (+ f_8_7 f_8_2 f_8_3 (- f_11_8)) (- 6.0)))  ;; node 8
(assert (= (+ f_9_6 (- f_10_9) f_9_2) (- 10.0)))  ;; node 9
(assert (= f_10_9 (- 3.0)))  ;; node 10
(assert (= (+ f_11_4 (- f_12_11) f_11_7 f_11_6 (- f_13_11) f_11_8 f_11_2) (- 10.0)))  ;; node 11
(assert (= f_12_11 (- 25.0)))  ;; node 12
(assert (= f_13_11 3.0))  ;; node 13
(assert (= f_14_7 34.0))  ;; node 14

;; Sanity: sum_i P_i = 0

(check-sat)
(get-model)
