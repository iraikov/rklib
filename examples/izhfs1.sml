(*
   Izhikevich Fast Spiking (FS) model.
*)


structure IzhikevichFS1 =
struct

open RungeKutta

type IzhState = (real * real)

(* scaler :: real -> IzhState -> IzhState *)
fun scaler (s,(v,u)) = let open Real in (s*v, s*u) end

(* summer :: IzhState -> IzhState -> IzhState *)
fun summer ((v1,u1), (v2,u2)) = let open Real in (v1 + v2, u1 + u2) end


(* Parameters of the model *)

val params =
    {
      Vpeak = 25.0,
      k = 1.0,
      Vt = ~55.0,
      Vr = ~40.0,
      Vb = ~55.0,
      Cm = 20.0,
      a = 0.2,
      b = 0.025,
      c = ~45.0,
      I = 100.0
    }

fun deriv {Vpeak,k,Vt,Vr,Vb,Cm,a,b,c,I}
	  (t,(v,u)) =
    let open Real
        fun UU (v) = (if (v < Vb) then 0.0 else (b * (v - Vb) * (v - Vb) * (v - Vb)))
    in
        (((k * (v - Vr) * (v - Vt)) - u + I) / Cm,
         (a * ((UU v) - u)))
    end

fun driver tstep
           {Vpeak,k,Vt,Vr,Vb,Cm,a,b,c,I}
	   t 
           (v,u) =
    let open Real
    in
        if ((v - Vpeak) >= 0.0) 
        then (t+tstep,(c, (u + (tstep * a * ((b * (c - (b * (v - Vb) * (v - Vb) * (v - Vb)) )) - u)))))
        else (t+tstep,(v,u))
    end
    


val initial = (~65.0,~1.625)

val tstep = 0.001

val rk4b: IzhState stepper1 = make_rk4b()
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


fun showst (t, (v, u)) = 
    String.concat [(showReal t), " ", (showReal v) , " ", (showReal u) ]


fun solver (tmax,stepper) (t,st) =
  let 
      val (tn,stn) = (driver tstep params t) (stepper tstep (t,st))
  in
      putStrLn (showst (t,st));
      if tn > tmax
      then (putStrLn "# All done!"; (tn,stn))
      else solver (tmax,stepper) (tn,stn)
  end


val (tn,_) = solver (1000.0,make_stepper params) (0.0,initial)

end
