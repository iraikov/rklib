
(*
 Module for testing Runge-Kutta integration of ODEs
 Based on code by Uwe Hollerbach <uh@alumni.caltech.edu>


 Solve the restricted 3-body problem, seek Arenstorf orbits

	x1'' = x1 + 2 x2' - mu'*(x1 + mu)/D1 - mu*(x1 - mu')/D2
	x2'' = x2 - 2 x1' - mu'*x2/D1 - mu*x2/D2
	D1 = ((x1 + mu)^2 + x2^2)^(3/2)
	D2 = ((x1 - mu')^2 + x2^2)^(3/2)
	mu' = 1 - mu

 Rewrite the second derivatives in terms of new variables v1 and v2 to
 get a system of first-order equations

 Represent the solution as a tuple (x1,v1,x2,v2), then we need three
 functions to manipulate these: a summer, a scaler, and a derivative
 operator: all very trivial functions.
*)

structure Arenstorf =
struct

open RungeKutta

type AState = (real * real * real * real)

(* scaler :: real -> AState -> AState *)
fun scaler (s,(x1,v1,x2,v2)) = let open Real in (s*x1, s*v1, s*x2, s*v2) end

(* summer :: AState -> AState -> AState *)
fun summer ((a1,a2,a3,a4),(b1,b2,b3,b4)) = let open Real in (a1+b1, a2+b2, a3+b3, a4+b4) end


(* Parameters of the model *)

val mu  = 0.012277471
val mu' = 1.0 - mu

(* from Hairer Norsett Wanner p 130: 3-loop orbit *)
val tmax3 = 17.0652165601579625588917206249
val start3 = (0.994, 0.0, 0.0, ~2.00158510637908252240537862224)

(* from Hairer Norsett Wanner p 186: 2-loop orbit *)
val tmax2 = 11.124340337266085134999734047
val start2 = (0.994, 0.0, 0.0, ~2.0317326295573368357302057924)

fun deriv (t,(x1,v1,x2,v2)) =
  let open Real

      val yp = x1 + mu
      val ym = x1 - mu'
      val d1 = Math.sqrt (yp*yp + x2*x2)
      val d1c = d1*d1*d1
      val d2 = Math.sqrt (ym*ym + x2*x2)
      val d2c = d2*d2*d2
      val f1 = x1 + 2.0*v2 - (mu'*yp/d1c) - (mu*ym/d2c)
      val f2 = x2 - 2.0*v1 - (mu'*x2/d1c) - (mu*x2/d2c)
  in (v1,f1,v2,f2) end

(* Step-size selector oracle *)

datatype ('a, 'b) either = Left of 'a
			 | Right of 'b

val tol = 1.0 * Real.Math.pow (10.0, ~9.0)

val lb = 0.5 * tol
val ub = 1.0 * tol

(* predictor :: real -> LVState ->  real real either *)
fun predictor (h,(x1,v1,x2,v2)) =
  let open Real
      val e = (abs x1) + (abs v1) + (abs x2) + (abs v2)
  in if e < lb 
        then Right (1.414*h)		(* step too small, accept but grow *)
        else if e < ub 
             then Right h		(* step just right *)
             else Left (0.5*h)		(* step too large, reject and shrink *)
  end



fun putStrLn str =
    (TextIO.output (TextIO.stdOut, str);
     TextIO.output (TextIO.stdOut, "\n"))

fun putStr str =
    (TextIO.output (TextIO.stdOut, str))

fun showReal n = Real.toString n

fun showpos (t,(x1,v1,x2,v2)) = String.concat [(showReal t) ^ " " ^ (showReal x1) ^ " " ^ (showReal x2)]
fun showst  (t,(x1,v1,x2,v2)) = String.concat [(showReal t), " ", (showReal x1) , " ", (showReal v1),
					       (showReal x2) , " ", (showReal v2)]


(* Numerical parameters: method we want to use, initial stepsize *)

val rkf78: AState stepper2  = make_rkf78() 
val stepper = rkf78 (scaler,summer,deriv)

val step = 100.0

fun solver tmax (h,t,st) =
  let open Real
      val hf = tmax - t
      val hs = if hf > h then h else hf
      val (tn, stn, etn) = stepper hs (t,st)

      fun shs (pre, h) = TextIO.output (TextIO.stdErr, (pre ^ "creasing step to " ^ (showReal h) ^ "\n"))
      fun shnil () = TextIO.output (TextIO.stdErr, "")
  in
      if hs < 1.0 andalso Real.== (hs, hf)
      then (putStrLn (showst (tn, stn)); putStrLn "# All done!")
      else
	  case predictor (h,etn) of
              Left bad => (shs ("de", bad); solver tmax (bad,t,st))
	    | Right good => ((if good > h
			      then shs ("in", good)
			      else shnil());
			     putStrLn (showst (t,st));
			     solver tmax (good,tn,stn))
  end


val _ = solver tmax2 (step,0.0,start2)
val _ = solver tmax3 (step,0.0,start3)

end


