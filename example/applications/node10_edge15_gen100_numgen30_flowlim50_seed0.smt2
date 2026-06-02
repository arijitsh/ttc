;; SMT-LRA encoding for power-flow toy model: node10_edge15_gen100_numgen30_flowlim50_seed0
;; Constants inlined; one flow variable per undirected line; f_j_i is -f_i_j.
;; Nodes: 10; Edges: 15

(set-logic QF_LRA)

;; ---------- Edge variables (canonical orientation i>j) ----------
;; Edge {6,5} with capacity 48
(declare-fun proj_x_6_5 () Bool)
(declare-fun f_6_5 () Real)
;; Edge {10,6} with capacity 50
(declare-fun proj_x_10_6 () Bool)
(declare-fun f_10_6 () Real)
;; Edge {10,7} with capacity 47
(declare-fun proj_x_10_7 () Bool)
(declare-fun f_10_7 () Real)
;; Edge {6,1} with capacity 47
(declare-fun proj_x_6_1 () Bool)
(declare-fun f_6_1 () Real)
;; Edge {10,2} with capacity 51
(declare-fun proj_x_10_2 () Bool)
(declare-fun f_10_2 () Real)
;; Edge {8,7} with capacity 49
(declare-fun proj_x_8_7 () Bool)
(declare-fun f_8_7 () Real)
;; Edge {3,1} with capacity 52
(declare-fun proj_x_3_1 () Bool)
(declare-fun f_3_1 () Real)
;; Edge {9,5} with capacity 51
(declare-fun proj_x_9_5 () Bool)
(declare-fun f_9_5 () Real)
;; Edge {4,2} with capacity 46
(declare-fun proj_x_4_2 () Bool)
(declare-fun f_4_2 () Real)
;; Edge {9,2} with capacity 49
(declare-fun proj_x_9_2 () Bool)
(declare-fun f_9_2 () Real)
;; Edge {7,4} with capacity 48
(declare-fun proj_x_7_4 () Bool)
(declare-fun f_7_4 () Real)
;; Edge {9,7} with capacity 49
(declare-fun proj_x_9_7 () Bool)
(declare-fun f_9_7 () Real)
;; Edge {5,4} with capacity 50
(declare-fun proj_x_5_4 () Bool)
(declare-fun f_5_4 () Real)
;; Edge {8,3} with capacity 52
(declare-fun proj_x_8_3 () Bool)
(declare-fun f_8_3 () Real)
;; Edge {8,4} with capacity 47
(declare-fun proj_x_8_4 () Bool)
(declare-fun f_8_4 () Real)

;; ---------- Line-outage and capacity limits ----------
(assert (=> (not proj_x_6_5) (= f_6_5 0.0)))
(assert (and (<= (- 48.0) f_6_5) (<= f_6_5 48.0)))
(assert (=> (not proj_x_10_6) (= f_10_6 0.0)))
(assert (and (<= (- 50.0) f_10_6) (<= f_10_6 50.0)))
(assert (=> (not proj_x_10_7) (= f_10_7 0.0)))
(assert (and (<= (- 47.0) f_10_7) (<= f_10_7 47.0)))
(assert (=> (not proj_x_6_1) (= f_6_1 0.0)))
(assert (and (<= (- 47.0) f_6_1) (<= f_6_1 47.0)))
(assert (=> (not proj_x_10_2) (= f_10_2 0.0)))
(assert (and (<= (- 51.0) f_10_2) (<= f_10_2 51.0)))
(assert (=> (not proj_x_8_7) (= f_8_7 0.0)))
(assert (and (<= (- 49.0) f_8_7) (<= f_8_7 49.0)))
(assert (=> (not proj_x_3_1) (= f_3_1 0.0)))
(assert (and (<= (- 52.0) f_3_1) (<= f_3_1 52.0)))
(assert (=> (not proj_x_9_5) (= f_9_5 0.0)))
(assert (and (<= (- 51.0) f_9_5) (<= f_9_5 51.0)))
(assert (=> (not proj_x_4_2) (= f_4_2 0.0)))
(assert (and (<= (- 46.0) f_4_2) (<= f_4_2 46.0)))
(assert (=> (not proj_x_9_2) (= f_9_2 0.0)))
(assert (and (<= (- 49.0) f_9_2) (<= f_9_2 49.0)))
(assert (=> (not proj_x_7_4) (= f_7_4 0.0)))
(assert (and (<= (- 48.0) f_7_4) (<= f_7_4 48.0)))
(assert (=> (not proj_x_9_7) (= f_9_7 0.0)))
(assert (and (<= (- 49.0) f_9_7) (<= f_9_7 49.0)))
(assert (=> (not proj_x_5_4) (= f_5_4 0.0)))
(assert (and (<= (- 50.0) f_5_4) (<= f_5_4 50.0)))
(assert (=> (not proj_x_8_3) (= f_8_3 0.0)))
(assert (and (<= (- 52.0) f_8_3) (<= f_8_3 52.0)))
(assert (=> (not proj_x_8_4) (= f_8_4 0.0)))
(assert (and (<= (- 47.0) f_8_4) (<= f_8_4 47.0)))

;; ---------- Nodal balance ----------
(assert (= (+ (- f_6_1) (- f_3_1)) 34.0))  ;; node 1
(assert (= (+ (- f_10_2) (- f_4_2) (- f_9_2)) (- 39.0)))  ;; node 2
(assert (= (+ f_3_1 (- f_8_3)) (- 7.0)))  ;; node 3
(assert (= (+ f_4_2 (- f_7_4) (- f_5_4) (- f_8_4)) (- 6.0)))  ;; node 4
(assert (= (+ (- f_6_5) (- f_9_5) f_5_4) (- 10.0)))  ;; node 5
(assert (= (+ f_6_5 (- f_10_6) f_6_1) (- 1.0)))  ;; node 6
(assert (= (+ (- f_10_7) (- f_8_7) f_7_4 (- f_9_7)) 32.0))  ;; node 7
(assert (= (+ f_8_7 f_8_3 f_8_4) (- 12.0)))  ;; node 8
(assert (= (+ f_9_5 f_9_2 f_9_7) (- 25.0)))  ;; node 9
(assert (= (+ f_10_6 f_10_7 f_10_2) 34.0))  ;; node 10

;; Sanity: sum_i P_i = 0

(check-sat)
(get-model)
