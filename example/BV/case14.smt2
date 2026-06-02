; Variables:
;   Xb[i]    : bus i enabled
;   theta[i] : voltage angle at bus i
;   Xg[j]    : generator j enabled
;   Pg[j]    : power output of generator j
;   projXl[k]    : line k enabled
;   F[k]     : power flow on line k

; Constraints:
;   generator limits and zero output when off
;   line flow equations and thermal limits
;   power balance at every bus

(set-logic QF_LRA)
; benchmark generated from python API
(set-info :status unknown)
(declare-fun theta1 () Real)
(declare-fun Pg0 () Real)
(declare-fun Xg0 () Bool)
(declare-fun Pg1 () Real)
(declare-fun Xg1 () Bool)
(declare-fun Pg2 () Real)
(declare-fun Xg2 () Bool)
(declare-fun Pg3 () Real)
(declare-fun Xg3 () Bool)
(declare-fun Pg4 () Real)
(declare-fun Xg4 () Bool)
(declare-fun theta2 () Real)
(declare-fun F0 () Real)
(declare-fun projXl0 () Bool)
(declare-fun theta5 () Real)
(declare-fun F1 () Real)
(declare-fun projXl1 () Bool)
(declare-fun theta3 () Real)
(declare-fun F2 () Real)
(declare-fun projXl2 () Bool)
(declare-fun theta4 () Real)
(declare-fun F3 () Real)
(declare-fun projXl3 () Bool)
(declare-fun F4 () Real)
(declare-fun projXl4 () Bool)
(declare-fun F5 () Real)
(declare-fun projXl5 () Bool)
(declare-fun F6 () Real)
(declare-fun projXl6 () Bool)
(declare-fun theta7 () Real)
(declare-fun F7 () Real)
(declare-fun projXl7 () Bool)
(declare-fun theta9 () Real)
(declare-fun F8 () Real)
(declare-fun projXl8 () Bool)
(declare-fun theta6 () Real)
(declare-fun F9 () Real)
(declare-fun projXl9 () Bool)
(declare-fun theta11 () Real)
(declare-fun F10 () Real)
(declare-fun projXl10 () Bool)
(declare-fun theta12 () Real)
(declare-fun F11 () Real)
(declare-fun projXl11 () Bool)
(declare-fun theta13 () Real)
(declare-fun F12 () Real)
(declare-fun projXl12 () Bool)
(declare-fun theta8 () Real)
(declare-fun F13 () Real)
(declare-fun projXl13 () Bool)
(declare-fun F14 () Real)
(declare-fun projXl14 () Bool)
(declare-fun theta10 () Real)
(declare-fun F15 () Real)
(declare-fun projXl15 () Bool)
(declare-fun theta14 () Real)
(declare-fun F16 () Real)
(declare-fun projXl16 () Bool)
(declare-fun F17 () Real)
(declare-fun projXl17 () Bool)
(declare-fun F18 () Real)
(declare-fun projXl18 () Bool)
(declare-fun F19 () Real)
(declare-fun projXl19 () Bool)
(declare-fun Xb1 () Bool)
(declare-fun Xb2 () Bool)
(declare-fun Xb3 () Bool)
(declare-fun Xb4 () Bool)
(declare-fun Xb5 () Bool)
(declare-fun Xb6 () Bool)
(declare-fun Xb7 () Bool)
(declare-fun Xb8 () Bool)
(declare-fun Xb9 () Bool)
(declare-fun Xb10 () Bool)
(declare-fun Xb11 () Bool)
(declare-fun Xb12 () Bool)
(declare-fun Xb13 () Bool)
(declare-fun Xb14 () Bool)
(assert
 (= theta1 0.0))
(assert
 (let (($x106 (and (>= Pg0 (/ 1662.0 5.0)) (<= Pg0 0.0))))
 (=> Xg0 $x106)))
(assert
 (let (($x110 (= Pg0 0.0)))
 (let (($x108 (not Xg0)))
 (=> $x108 $x110))))
(assert
 (let (($x116 (and (>= Pg1 140.0) (<= Pg1 0.0))))
 (=> Xg1 $x116)))
(assert
 (let (($x120 (= Pg1 0.0)))
 (let (($x118 (not Xg1)))
 (=> $x118 $x120))))
(assert
 (let (($x126 (and (>= Pg2 100.0) (<= Pg2 0.0))))
 (=> Xg2 $x126)))
(assert
 (let (($x130 (= Pg2 0.0)))
 (let (($x128 (not Xg2)))
 (=> $x128 $x130))))
(assert
 (let (($x135 (and (>= Pg3 100.0) (<= Pg3 0.0))))
 (=> Xg3 $x135)))
(assert
 (let (($x139 (= Pg3 0.0)))
 (let (($x137 (not Xg3)))
 (=> $x137 $x139))))
(assert
 (let (($x144 (and (>= Pg4 100.0) (<= Pg4 0.0))))
 (=> Xg4 $x144)))
(assert
 (let (($x148 (= Pg4 0.0)))
 (let (($x146 (not Xg4)))
 (=> $x146 $x148))))
(assert
 (= F0 (* (/ 5917.0 100000.0) (- theta1 theta2))))
(assert
 (let (($x171 (and (<= F0 9900.0) (>= F0 (- 9900.0)))))
 (=> projXl0 $x171)))
(assert
 (let (($x175 (= F0 0.0)))
 (let (($x173 (not projXl0)))
 (=> $x173 $x175))))
(assert
 (= F1 (* (/ 697.0 3125.0) (- theta1 theta5))))
(assert
 (let (($x195 (and (<= F1 9900.0) (>= F1 (- 9900.0)))))
 (=> projXl1 $x195)))
(assert
 (let (($x199 (= F1 0.0)))
 (let (($x197 (not projXl1)))
 (=> $x197 $x199))))
(assert
 (= F2 (* (/ 19797.0 100000.0) (- theta2 theta3))))
(assert
 (let (($x219 (and (<= F2 9900.0) (>= F2 (- 9900.0)))))
 (=> projXl2 $x219)))
(assert
 (let (($x223 (= F2 0.0)))
 (let (($x221 (not projXl2)))
 (=> $x221 $x223))))
(assert
 (= F3 (* (/ 551.0 3125.0) (- theta2 theta4))))
(assert
 (let (($x243 (and (<= F3 9900.0) (>= F3 (- 9900.0)))))
 (=> projXl3 $x243)))
(assert
 (let (($x247 (= F3 0.0)))
 (let (($x245 (not projXl3)))
 (=> $x245 $x247))))
(assert
 (= F4 (* (/ 4347.0 25000.0) (- theta2 theta5))))
(assert
 (let (($x266 (and (<= F4 9900.0) (>= F4 (- 9900.0)))))
 (=> projXl4 $x266)))
(assert
 (let (($x270 (= F4 0.0)))
 (let (($x268 (not projXl4)))
 (=> $x268 $x270))))
(assert
 (= F5 (* (/ 17103.0 100000.0) (- theta3 theta4))))
(assert
 (let (($x289 (and (<= F5 9900.0) (>= F5 (- 9900.0)))))
 (=> projXl5 $x289)))
(assert
 (let (($x293 (= F5 0.0)))
 (let (($x291 (not projXl5)))
 (=> $x291 $x293))))
(assert
 (= F6 (* (/ 4211.0 100000.0) (- theta4 theta5))))
(assert
 (let (($x312 (and (<= F6 9900.0) (>= F6 (- 9900.0)))))
 (=> projXl6 $x312)))
(assert
 (let (($x316 (= F6 0.0)))
 (let (($x314 (not projXl6)))
 (=> $x314 $x316))))
(assert
 (= F7 (* (/ 1307.0 6250.0) (- theta4 theta7))))
(assert
 (let (($x336 (and (<= F7 9900.0) (>= F7 (- 9900.0)))))
 (=> projXl7 $x336)))
(assert
 (let (($x340 (= F7 0.0)))
 (let (($x338 (not projXl7)))
 (=> $x338 $x340))))
(assert
 (= F8 (* (/ 27809.0 50000.0) (- theta4 theta9))))
(assert
 (let (($x360 (and (<= F8 9900.0) (>= F8 (- 9900.0)))))
 (=> projXl8 $x360)))
(assert
 (let (($x364 (= F8 0.0)))
 (let (($x362 (not projXl8)))
 (=> $x362 $x364))))
(assert
 (= F9 (* (/ 12601.0 50000.0) (- theta5 theta6))))
(assert
 (let (($x384 (and (<= F9 9900.0) (>= F9 (- 9900.0)))))
 (=> projXl9 $x384)))
(assert
 (let (($x388 (= F9 0.0)))
 (let (($x386 (not projXl9)))
 (=> $x386 $x388))))
(assert
 (= F10 (* (/ 1989.0 10000.0) (- theta6 theta11))))
(assert
 (let (($x408 (and (<= F10 9900.0) (>= F10 (- 9900.0)))))
 (=> projXl10 $x408)))
(assert
 (let (($x412 (= F10 0.0)))
 (let (($x410 (not projXl10)))
 (=> $x410 $x412))))
(assert
 (= F11 (* (/ 25581.0 100000.0) (- theta6 theta12))))
(assert
 (let (($x432 (and (<= F11 9900.0) (>= F11 (- 9900.0)))))
 (=> projXl11 $x432)))
(assert
 (let (($x436 (= F11 0.0)))
 (let (($x434 (not projXl11)))
 (=> $x434 $x436))))
(assert
 (= F12 (* (/ 13027.0 100000.0) (- theta6 theta13))))
(assert
 (let (($x456 (and (<= F12 9900.0) (>= F12 (- 9900.0)))))
 (=> projXl12 $x456)))
(assert
 (let (($x460 (= F12 0.0)))
 (let (($x458 (not projXl12)))
 (=> $x458 $x460))))
(assert
 (= F13 (* (/ 3523.0 20000.0) (- theta7 theta8))))
(assert
 (let (($x480 (and (<= F13 9900.0) (>= F13 (- 9900.0)))))
 (=> projXl13 $x480)))
(assert
 (let (($x484 (= F13 0.0)))
 (let (($x482 (not projXl13)))
 (=> $x482 $x484))))
(assert
 (= F14 (* (/ 11001.0 100000.0) (- theta7 theta9))))
(assert
 (let (($x503 (and (<= F14 9900.0) (>= F14 (- 9900.0)))))
 (=> projXl14 $x503)))
(assert
 (let (($x507 (= F14 0.0)))
 (let (($x505 (not projXl14)))
 (=> $x505 $x507))))
(assert
 (= F15 (* (/ 169.0 2000.0) (- theta9 theta10))))
(assert
 (let (($x527 (and (<= F15 9900.0) (>= F15 (- 9900.0)))))
 (=> projXl15 $x527)))
(assert
 (let (($x531 (= F15 0.0)))
 (let (($x529 (not projXl15)))
 (=> $x529 $x531))))
(assert
 (= F16 (* (/ 13519.0 50000.0) (- theta9 theta14))))
(assert
 (let (($x551 (and (<= F16 9900.0) (>= F16 (- 9900.0)))))
 (=> projXl16 $x551)))
(assert
 (let (($x555 (= F16 0.0)))
 (let (($x553 (not projXl16)))
 (=> $x553 $x555))))
(assert
 (= F17 (* (/ 19207.0 100000.0) (- theta10 theta11))))
(assert
 (let (($x574 (and (<= F17 9900.0) (>= F17 (- 9900.0)))))
 (=> projXl17 $x574)))
(assert
 (let (($x578 (= F17 0.0)))
 (let (($x576 (not projXl17)))
 (=> $x576 $x578))))
(assert
 (= F18 (* (/ 4997.0 25000.0) (- theta12 theta13))))
(assert
 (let (($x597 (and (<= F18 9900.0) (>= F18 (- 9900.0)))))
 (=> projXl18 $x597)))
(assert
 (let (($x601 (= F18 0.0)))
 (let (($x599 (not projXl18)))
 (=> $x599 $x601))))
(assert
 (= F19 (* (/ 17401.0 50000.0) (- theta13 theta14))))
(assert
 (let (($x620 (and (<= F19 9900.0) (>= F19 (- 9900.0)))))
 (=> projXl19 $x620)))
(assert
 (let (($x624 (= F19 0.0)))
 (let (($x622 (not projXl19)))
 (=> $x622 $x624))))
(assert
 (=> Xb1 (= (- Pg0 0.0) (+ F0 F1))))
(assert
 (=> Xb2 (= (- Pg1 (/ 217.0 10.0)) (+ (- F0) F2 F3 F4))))
(assert
 (=> Xb3 (= (- Pg2 (/ 471.0 5.0)) (+ (- F2) F5))))
(assert
 (=> Xb4 (= (- 0.0 (/ 239.0 5.0)) (+ (- F3) (- F5) F6 F7 F8))))
(assert
 (=> Xb5 (= (- 0.0 (/ 38.0 5.0)) (+ (- F1) (- F4) (- F6) F9))))
(assert
 (=> Xb6 (= (- Pg3 (/ 56.0 5.0)) (+ (- F9) F10 F11 F12))))
(assert
 (=> Xb7 (= (- 0.0 0.0) (+ (- F7) F13 F14))))
(assert
 (=> Xb8 (= (- Pg4 0.0) (- F13))))
(assert
 (=> Xb9 (= (- 0.0 (/ 59.0 2.0)) (+ (- F8) (- F14) F15 F16))))
(assert
 (=> Xb10 (= (- 0.0 9.0) (+ (- F15) F17))))
(assert
 (=> Xb11 (= (- 0.0 (/ 7.0 2.0)) (+ (- F10) (- F17)))))
(assert
 (=> Xb12 (= (- 0.0 (/ 61.0 10.0)) (+ (- F11) F18))))
(assert
 (=> Xb13 (= (- 0.0 (/ 27.0 2.0)) (+ (- F12) (- F18) F19))))
(assert
 (=> Xb14 (= (- 0.0 (/ 149.0 10.0)) (+ (- F16) (- F19)))))
(check-sat)
