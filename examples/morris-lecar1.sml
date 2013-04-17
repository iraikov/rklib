(*
   Morris-Lecar reduced model.

   for type I dynamics, use v3=12,v4=17,phi=.06666667 
   for type II dynamics, use v3=2,v4=30,phi=.04
*)


structure MorrisLecar1 =
struct

open RungeKutta

type MLState = (real * real)

(* scaler :: real -> MLState -> MLState *)
fun scaler (s,(y1,y2)) = let open Real in (s*y1, s*y2) end

(* summer :: MLState -> MLState -> MLState *)
fun summer ((y11,y12), (y21,y22)) = let open Real in (y11 + y21, y12 + y22) end


(* Parameters of the model *)

val type1 =
    {Istim = 20.0,
     vk    = ~84.0,
     vl    = ~60.0,
     vca   = 120.0,
     gk    = 2.5,
     gl    = 0.5,
     gca   = 1.0,
     c     = 2.0,
     v1    = ~1.2,
     v2    = 18.0,
     v3    = 12.0,
     v4    = 17.0,
     phi   = 0.0667}
		    
val type2 =
    {Istim = 20.0,
     vk    = ~84.0,
     vl    = ~60.0,
     vca   = 120.0,
     gk    = 2.5,
     gl    = 0.5,
     gca   = 1.0,
     c     = 2.0,
     v1    = ~1.2,
     v2    = 18.0,
     v3    = 2.0,
     v4    = 30.0,
     phi   = 0.04}
		    

fun deriv {Istim, vk,vl,vca,gk,gl,gca,c,v1,v2,v3,v4,phi}
	  (t,(v,w)) =
    let open Real
	fun minf (v) = 0.5 * (1.0 + (Math.tanh ((v - v1) / v2))) 
	fun winf (v) = 0.5 * (1.0 + (Math.tanh ((v - v3) / v4))) 
	fun lamw (v) = phi * (Math.cosh ((v / v3) / (2.0 * v4)))
		    
	val ica = gca * ((minf v) *  (vca - v))
	val ik  = gk  * (w * (vk - v))
    in
	( (Istim + (gl * (vl - v)) + ica + ik) / c,
	  (lamw v) * ((winf v) - w) )
    end

val initial = (~60.899,0.0149)

val step = 0.001

val rk4b: MLState stepper1 = make_rk4b()
fun make_stepper (params) = rk4b (scaler,summer,deriv params)


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

fun solver (tmax,stepper) (t,st) =
  let val (tn,stn) = stepper step (t,st)
  in putStrLn (showst (t,st));
     if t > tmax
     then (putStrLn "# All done!"; (tn,stn))
     else solver (tmax,stepper) (tn,stn)
  end


val (tn,_) = solver (100.0,make_stepper type1) (0.0,initial)

end
