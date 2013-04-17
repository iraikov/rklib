
(*
 Module for testing Runge-Kutta integration of ODEs
 Based on code by , Uwe Hollerbach <uh@alumni.caltech.edu>.


 Solve the Lotka-Volterra equation

 dy1                   dy2
 --- = y1*(a - b*y2)   --- = -y2*(c - d*y1)
  dt                    dt

 Represent the solution as a tuple (y1,y2), then we need three functions
 to manipulate these: a summer, a scaler, and a derivative operator: all
 very trivial functions.

 Y1 is the prey population, Y2 the predator
*)

structure LotkaVolterra2 =
struct

open RungeKutta

type LVState = (real * real)

(* scaler :: real -> LVState -> LVState *)
fun scaler (s,(y1,y2)) = let open Real in (s*y1, s*y2) end

(* summer :: LVState -> LVState -> LVState *)
fun summer ((y11,y12), (y21,y22)) = let open Real in (y11 + y21, y12 + y22) end

(* Parameters of the model *)

val a = 3.0    (* prey population growth parameter *)
val b = 2.0    (* prey decay due to presence of predator *)
val c = 3.0    (* predator starvation parameter *)
val d = 2.0    (* predator growth due to presence of prey *)

(* deriv :: real -> LVState -> LVState *)
fun deriv (t,(y1,y2)) = let open Real in (y1*(a - b*y2), ~y2*(c - d*y1)) end

val rkdp: LVState stepper2 = make_rkdp()
val stepper = rkdp (scaler,summer,deriv)

datatype ('a, 'b) either = Left of 'a
			 | Right of 'b

val lb = 0.5 * Real.Math.pow (10.0, ~8.0)
val ub = 1.0 * Real.Math.pow (10.0, ~8.0)

(* predictor :: real -> LVState ->  real real either *)
fun predictor (h,(y1,y2)) =
  let open Real
      val e = (abs y1) + (abs y2)
  in if e < lb 
        then Right (1.414*h)		(* step too small, accept but grow *)
        else if e < ub 
             then Right h		(* step just right *)
             else Left (0.5*h)		(* step too large, reject and shrink *)
  end

(* initial step guaranteed to fail *)
val step = 100.0

fun putStrLn str =
    (TextIO.output (TextIO.stdOut, str);
     TextIO.output (TextIO.stdOut, "\n"))

fun putStr str =
    (TextIO.output (TextIO.stdOut, str))

fun showReal n = Real.toString n

fun showst (t, (y1, y2)) = String.concat [(showReal t), " ", (showReal y1) , " ", (showReal y2) ]

fun solver (h,t,st) =
  let open Real
      val hf = 40.0 - t
      val hs = if hf > h then h else hf
      val (tn, stn, etn) = stepper hs (t,st)

      fun shs (pre, h) = TextIO.output (TextIO.stdErr, (pre ^ "creasing step to " ^ (showReal h) ^ "\n"))
      fun shnil () = TextIO.output (TextIO.stdErr, "")
  in
      if hs < 1.0 andalso Real.== (hs, hf)
      then (putStrLn (showst (tn, stn)); putStrLn "# All done!")
      else
	  case predictor (h,etn) of
              Left bad => (shs ("de", bad); solver (bad,t,st))
	    | Right good => ((if good > h
			      then shs ("in", good)
			      else shnil());
			     putStrLn (showst (t,st));
			     solver (good,tn,stn))
  end


val _ = solver (step,0.0,(2.0,0.3))

end

