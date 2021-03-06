(*
   Morris-Lecar reduced model.

   for type I dynamics, use v3=12,v4=17,phi=.06666667 
   for type II dynamics, use v3=2,v4=30,phi=.04
*)


structure MorrisLecar2 =
struct

open RungeKutta

type MLState = (real * real)

(* scaler :: real -> MLState -> MLState *)
fun scaler (s,(y1,y2)) = let open Real in (s*y1, s*y2) end

(* summer :: MLState -> MLState -> MLState *)
fun summer ((y11,y12), (y21,y22)) = let open Real in (y11 + y21, y12 + y22) end


(* Parameters of the model *)

val params =
    {Istim = 50.0,
     c     = 20.0,
     vk    = ~70.0,
     vl    = ~50.0,
     vca   = 100.0,
     gk    = 8.0,
     gl    = 2.0,
     gca   = 4.0,
     v1    = ~1.0,
     v2    = 15.0,
     v3    = 10.0,
     v4    = 14.0,
     phi   = 0.0667}
		    
		    
fun deriv {Istim, vk,vl,vca,gk,gl,gca,c,v1,v2,v3,v4,phi}
	  (t,(v,w)) =
    let open Real
	fun minf (v) = 0.5 * (1.0 + (Math.tanh ((v - v1) / v2))) 
	fun winf (v) = 0.5 * (1.0 + (Math.tanh ((v - v3) / v4))) 
	fun lamw (v) = phi * (Math.cosh ((v - v3) / (2.0 * v4)))
		    
	val ica = gca * ((minf v) *  (vca - v))
	val ik  = gk  * (w * (vk - v))
    in
	( (Istim + (gl * (vl - v)) + ica + ik) / c,
	  (lamw v) * ((winf v) - w) )
    end

val initial: MLState = (~60.899,0.0149)

val step = 0.01

val rkdp: MLState stepper2 = make_rkdp()
val rkf45: MLState stepper2 = make_rkf45()
val rkf78: MLState stepper2 = make_rkf78()

fun make_stepper (params) = rkdp (scaler,summer,deriv params)

datatype ('a, 'b) either = Left of 'a | Right of 'b

val tol = Real.Math.pow (10.0, ~6.0)
val lb = 0.5 * tol
val ub = 1.0 * tol

(* predictor :: real -> real * MLState ->  real real either *)
fun predictor tol (h,(v,w)) =
  let open Real
      val e = (abs v) + (abs w)
  in 
      if e < lb 
      then Right (1.414*h)		(* step too small, accept but grow *)
      else (if e < ub 
            then Right h		(* step just right *)
            else Left (0.5*h))		(* step too large, reject and shrink *)
  end

fun putStrLn str =
    (TextIO.output (TextIO.stdOut, str);
     TextIO.output (TextIO.stdOut, "\n"))

fun putStr str =
    (TextIO.output (TextIO.stdOut, str))

fun showReal n = 
    let open Real 
	open StringCvt
    in
	(if n < 0.0 then "-" else "") ^ (Real.fmt (FIX (SOME 12)) (Real.abs n))
    end

fun showst (t, (y1, y2)) = String.concat [(showReal t), " ", (showReal y1) , " ", (showReal y2) ]

fun solver (tmax,stepper) (h,t,st) =
  let open Real

      val hf = tmax - t
      val hs = if hf > h then h else hf
      val (stn, etn) = stepper hs (t,st)

      fun shs (pre, h) = TextIO.output (TextIO.stdErr, ("# " ^ pre ^ "creasing step to " ^ (showReal h) ^ "\n"))

  in
      if hs < 1.0 andalso Real.== (hs, hf)
      then (putStrLn (showst (tmax, stn)); putStrLn "# All done!"; (tmax, stn))
      else
	  case predictor tol (h,etn) of
              Left bad => (shs ("de", bad); solver (tmax,stepper) (bad,t,st))
	    | Right good => 
              let 
                  val tn = t + hs
              in
                  ((if good > h
		    then shs ("in", good)
		    else ());
	           putStrLn (showst (t,st));
	           solver (tmax,stepper) (good,tn,stn))
              end
  end

val tmax = 1000.0

val stepper1 = make_stepper params

val (tn,_) = solver (tmax,stepper1) (step,0.0,initial)


end
