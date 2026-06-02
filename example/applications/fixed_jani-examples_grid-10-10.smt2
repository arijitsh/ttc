(set-logic QF_LIA)
(declare-fun target_10 () Bool)
(declare-fun target_9 () Bool)
(declare-fun target_8 () Bool)
(declare-fun target_7 () Bool)
(declare-fun target_6 () Bool)
(declare-fun target_5 () Bool)
(declare-fun target_4 () Bool)
(declare-fun target_3 () Bool)
(declare-fun target_2 () Bool)
(declare-fun target_1 () Bool)
(declare-fun target_0 () Bool)
(declare-fun x_10 () Int)
(declare-fun x_9 () Int)
(declare-fun x_8 () Int)
(declare-fun x_7 () Int)
(declare-fun x_6 () Int)
(declare-fun x_5 () Int)
(declare-fun x_4 () Int)
(declare-fun x_3 () Int)
(declare-fun x_2 () Int)
(declare-fun x_1 () Int)
(declare-fun x_0 () Int)
(declare-fun proj-0-1_9 () Bool)
(declare-fun proj-0-0_9 () Bool)
(declare-fun y_9 () Int)
(declare-fun y_10 () Int)
(declare-fun proj-0-1_8 () Bool)
(declare-fun proj-0-0_8 () Bool)
(declare-fun y_8 () Int)
(declare-fun proj-0-1_7 () Bool)
(declare-fun proj-0-0_7 () Bool)
(declare-fun y_7 () Int)
(declare-fun proj-0-1_6 () Bool)
(declare-fun proj-0-0_6 () Bool)
(declare-fun y_6 () Int)
(declare-fun proj-0-1_5 () Bool)
(declare-fun proj-0-0_5 () Bool)
(declare-fun y_5 () Int)
(declare-fun proj-0-1_4 () Bool)
(declare-fun proj-0-0_4 () Bool)
(declare-fun y_4 () Int)
(declare-fun proj-0-1_3 () Bool)
(declare-fun proj-0-0_3 () Bool)
(declare-fun y_3 () Int)
(declare-fun proj-0-1_2 () Bool)
(declare-fun proj-0-0_2 () Bool)
(declare-fun y_2 () Int)
(declare-fun proj-0-1_1 () Bool)
(declare-fun proj-0-0_1 () Bool)
(declare-fun y_1 () Int)
(declare-fun proj-0-1_0 () Bool)
(declare-fun proj-0-0_0 () Bool)
(declare-fun y_0 () Int)
(assert (let ((a!1 (=> (and (and (< x_0 10) (< y_0 10)) proj-0-0_0)
               (and (= x_1 (+ x_0 1)) (= y_1 y_0))))
      (a!2 (=> (and (and (< x_0 10) (< y_0 10)) proj-0-1_0)
               (and (= y_1 (+ y_0 1)) (= x_1 x_0))))
      (a!3 (=> (and (and (< x_1 10) (< y_1 10)) proj-0-0_1)
               (and (= x_2 (+ x_1 1)) (= y_2 y_1))))
      (a!4 (=> (and (and (< x_1 10) (< y_1 10)) proj-0-1_1)
               (and (= y_2 (+ y_1 1)) (= x_2 x_1))))
      (a!5 (=> (and (and (< x_2 10) (< y_2 10)) proj-0-0_2)
               (and (= x_3 (+ x_2 1)) (= y_3 y_2))))
      (a!6 (=> (and (and (< x_2 10) (< y_2 10)) proj-0-1_2)
               (and (= y_3 (+ y_2 1)) (= x_3 x_2))))
      (a!7 (=> (and (and (< x_3 10) (< y_3 10)) proj-0-0_3)
               (and (= x_4 (+ x_3 1)) (= y_4 y_3))))
      (a!8 (=> (and (and (< x_3 10) (< y_3 10)) proj-0-1_3)
               (and (= y_4 (+ y_3 1)) (= x_4 x_3))))
      (a!9 (=> (and (and (< x_4 10) (< y_4 10)) proj-0-0_4)
               (and (= x_5 (+ x_4 1)) (= y_5 y_4))))
      (a!10 (=> (and (and (< x_4 10) (< y_4 10)) proj-0-1_4)
                (and (= y_5 (+ y_4 1)) (= x_5 x_4))))
      (a!11 (=> (and (and (< x_5 10) (< y_5 10)) proj-0-0_5)
                (and (= x_6 (+ x_5 1)) (= y_6 y_5))))
      (a!12 (=> (and (and (< x_5 10) (< y_5 10)) proj-0-1_5)
                (and (= y_6 (+ y_5 1)) (= x_6 x_5))))
      (a!13 (=> (and (and (< x_6 10) (< y_6 10)) proj-0-0_6)
                (and (= x_7 (+ x_6 1)) (= y_7 y_6))))
      (a!14 (=> (and (and (< x_6 10) (< y_6 10)) proj-0-1_6)
                (and (= y_7 (+ y_6 1)) (= x_7 x_6))))
      (a!15 (=> (and (and (< x_7 10) (< y_7 10)) proj-0-0_7)
                (and (= x_8 (+ x_7 1)) (= y_8 y_7))))
      (a!16 (=> (and (and (< x_7 10) (< y_7 10)) proj-0-1_7)
                (and (= y_8 (+ y_7 1)) (= x_8 x_7))))
      (a!17 (=> (and (and (< x_8 10) (< y_8 10)) proj-0-0_8)
                (and (= x_9 (+ x_8 1)) (= y_9 y_8))))
      (a!18 (=> (and (and (< x_8 10) (< y_8 10)) proj-0-1_8)
                (and (= y_9 (+ y_8 1)) (= x_9 x_8))))
      (a!19 (=> (and (and (< x_9 10) (< y_9 10)) proj-0-0_9)
                (and (= x_10 (+ x_9 1)) (= y_10 y_9))))
      (a!20 (=> (and (and (< x_9 10) (< y_9 10)) proj-0-1_9)
                (and (= y_10 (+ y_9 1)) (= x_10 x_9)))))
  (and (= 0 x_0)
       (= 0 y_0)
       a!1
       a!2

       (or (not proj-0-0_0) (not proj-0-1_0))
       (or proj-0-0_0 proj-0-1_0)
       a!3
       a!4

       (or (not proj-0-0_1) (not proj-0-1_1))
       (or proj-0-0_1 proj-0-1_1)
       a!5
       a!6

       (or (not proj-0-0_2) (not proj-0-1_2))
       (or proj-0-0_2 proj-0-1_2)
       a!7
       a!8

       (or (not proj-0-0_3) (not proj-0-1_3))
       (or proj-0-0_3 proj-0-1_3)
       a!9
       a!10

       (or (not proj-0-0_4) (not proj-0-1_4))
       (or proj-0-0_4 proj-0-1_4)
       a!11
       a!12

       (or (not proj-0-0_5) (not proj-0-1_5))
       (or proj-0-0_5 proj-0-1_5)
       a!13
       a!14

       (or (not proj-0-0_6) (not proj-0-1_6))
       (or proj-0-0_6 proj-0-1_6)
       a!15
       a!16

       (or (not proj-0-0_7) (not proj-0-1_7))
       (or proj-0-0_7 proj-0-1_7)
       a!17
       a!18

       (or (not proj-0-0_8) (not proj-0-1_8))
       (or proj-0-0_8 proj-0-1_8)
       a!19
       a!20

       (or (not proj-0-0_9) (not proj-0-1_9))
       (or proj-0-0_9 proj-0-1_9)
       (= target_0 (= x_0 10))
       (= target_1 (= x_1 10))
       (= target_2 (= x_2 10))
       (= target_3 (= x_3 10))
       (= target_4 (= x_4 10))
       (= target_5 (= x_5 10))
       (= target_6 (= x_6 10))
       (= target_7 (= x_7 10))
       (= target_8 (= x_8 10))
       (= target_9 (= x_9 10))
       (= target_10 (= x_10 10))
       (or target_0
           target_1
           target_2
           target_3
           target_4
           target_5
           target_6
           target_7
           target_8
           target_9
           target_10))))
(check-sat)