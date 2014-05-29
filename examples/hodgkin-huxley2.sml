
(*
 This model is used to calculate the membrane potential assuming
 some initial state. The calculation is based on sodium ion flow,
 potassium ion flow and leakage ion flow. (Hodgkin, A. L. and Huxley,
 A. F. (1952) "A Quantitative Description of Membrane Current and its
 Application to Conduction and Excitation in Nerve" Journal of
 Physiology 117: 500-544)
*)


structure HodgkinHuxley2 =
struct

open RungeKutta

type HHState = (real * real * real * real)

(* scaler :: real -> HHState -> HHState *)
fun scaler (s,(v,m,h,n)) = let open Real in (s*v, s*m, s*h, s*n) end

(* summer :: HHState -> HHState -> HHState *)
fun summer ((a1,a2,a3,a4),(b1,b2,b3,b4)) = let open Real in (a1+b1, a2+b2, a3+b3, a4+b4) end

(* Parameters of the model *)

val hh52 =
    {I_stim    = 20.0,
     C_m      = 1.0,
     E_Na     = 50.0,
     E_K      = ~77.0,
     E_L      = ~54.4,
     gbar_Na  = 120.0,
     gbar_K   = 36.0,
     g_L      = 0.3}

		    
fun deriv {I_stim,C_m,E_Na,E_K,E_L,gbar_Na,gbar_K,g_L}
	  (t,(v,m,h,n)) =
    let open Real
	val neg = ~
	val exp = Math.exp

	fun amf (v) =  0.1  *  ((v + 40.0) /  (1.0 - (exp ((neg (v + 40.0)) / 10.0))))
	fun bmf (v) =  4.0  *  (exp ((neg (v + 65.0)) / 18.0))
	fun ahf (v) =  0.07 *  (exp ((neg (v + 65.0)) / 20.0))
	fun bhf (v) =  1.0  /  (1.0 + (exp ((neg (v + 35.0)) / 10.0)))
	fun anf (v) =  0.01 * ((v + 55.0) / (1.0 - (exp ((neg (v + 55.0)) / 10.0))))
	fun bnf (v) =  0.125 *  (exp ((neg (v + 65.0)) / 80.0))
							     
	(* transition rates at current step *)
	val am = (amf v)
	val an = (anf v)
	val ah = (ahf v)
	val bm = (bmf v)
	val bn = (bnf v)
	val bh = (bhf v)
		 
	val g_Na = gbar_Na * (h * (m*m*m))
	val g_K  = gbar_K  * (n*n*n*n)

	(* currents *)
	val I_Na   = (v - E_Na) * g_Na
	val I_K    = (v - E_K)  * g_K
	val I_L    = g_L * (v - E_L)
		  
    in
	( (I_stim - I_L - I_Na - I_K) / C_m,
	  (am * (1.0 - m)) - (bm * m),
	  (ah * (1.0 - h)) - (bh * h),
	  (an * (1.0 - n)) - (bn * n) )
    end
	
val initial = (~65.0, 0.052, 0.596, 0.317)
		  


val step = 1.0

val rkdp: HHState stepper2 = make_rkdp()
fun make_stepper (params) = rkdp (scaler,summer,deriv params)

datatype ('a, 'b) either = Left of 'a
			 | Right of 'b

val tol = Real.Math.pow (10.0, ~12.0)
val lb = 0.5 * tol
val ub = 1.0 * tol

(* predictor :: real -> HHState ->  real real either *)
fun predictor tol (step,(v,m,h,n)) =
  let open Real
      val e = (abs v) + (abs m) + (abs h) + (abs n)
  in if e < lb 
        then Right (1.414*step)		(* step too small, accept but grow *)
        else if e < ub 
             then Right step		(* step just right *)
             else Left (0.5*step)	(* step too large, reject and shrink *)
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

fun showst  (t,(v,m,h,n)) = String.concat [(showReal t), " ", (showReal v) , " ", (showReal m),
					   (showReal h) , " ", (showReal n)]

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


val (tn,_) = solver (1000.0,make_stepper hh52) (step,0.0,initial)

end
