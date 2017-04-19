(*
   Leaky integrate-and-fire neuron with adaptive integration.
*)


structure LeakyIF =
struct

open RungeKutta

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


fun showst (t, V) = 
    String.concat [(showReal t), " ", (showReal V)]

fun scaler (s,v) = let open Real in s*v end

fun summer (v1,v2) = let open Real in v1 + v2 end


(* Parameters of the model *)

val params =
    {
      theta = 30.0,
      tau = 20.0,
      tau_rp = 5.0,
      Vreset = 0.0,
      R = 1.0,
      I = 50.0
    }


fun deriv {theta,tau,tau_rp,Vreset,R,I}
	  (t,V) =
    let
        open Real
    in
        ( ~V + I*R ) / tau
    end

fun make_event (params as {theta,tau,tau_rp,Vreset,R,I})  =
    let open Real 
    in
        (fn (V) => (V - theta),
         fn (V,tstep) => 
            let
                val V' = Vreset
            in
                V'
            end)
    end

datatype ('a, 'b) either = Left of 'a | Right of 'b

val tol = Real.Math.pow (10.0, ~6.0)
val lb = 0.5 * tol
val ub = 1.0 * tol

(* predictor :: real -> real * MLState ->  real real either *)
fun predictor tol (h,V) =
  let open Real
      val e = (abs V)
      val _ = putStrLn ("# predictor: e = " ^ (showReal e))

  in 
      if e < lb 
      then Right (1.414*h)		(* step too small, accept but grow *)
      else (if e < ub 
            then Right h		(* step just right *)
            else Left (0.5*h))		(* step too large, reject and shrink *)
  end



datatype 'a result = Next of 'a | Root of 'a
fun solver (stepper,evtest,hinterp) =
  let
      fun step (t,v,tstep) =
        let open Real
            val (v',ev,ks) = stepper tstep (t,v)
        in
            case predictor tol (tstep,ev) of
                Right tstep' => 
                (if (evtest (v') >= 0.0)
                 then (let
                          val tn      = t+tstep
                          val finterp = hinterp (tstep, ks, t, v)
                          val theta   = RootFind.brent tol (evtest o finterp) 0.0 1.0
                          val v'     = finterp theta
                      in
                          (putStrLn ("# event: theta = " ^ (showReal theta) ^ " tstep = " ^ (showReal tstep));
                           putStrLn ("# event: finterp theta = " ^ (showReal v'));
                           Root (t+tstep*theta,v',tstep'))
                      end)
                 else Next (t+tstep,v',tstep'))
              | Left tstep'  => 
                step (t,v,tstep')
        end
  in
      step
  end
    


val initial = 0.0

val tstep = 0.5

val cerkdp: real stepper3 = make_cerkdp()
fun make_stepper (params) = cerkdp (scaler,summer,deriv params)
val hinterp: real hinterp = make_interp_cerkdp (scaler,summer)

fun driver (tmax,stepper,(evtest,evhandle),hinterp) =
  let val solver' = solver (stepper,evtest,hinterp)
      fun recur (t,st,tstep) =
        let open Real in
            case solver' (t,st,tstep) of
                Next (tn,stn,tstep') =>
                (putStrLn (showst (t,st));
                 if tn > tmax
                 then (putStrLn "# All done!"; (tn,stn))
                 else recur (tn,stn,tstep'))
              | Root (tn,stn,tstep') =>
                (putStrLn (showst (t,st));
                 putStrLn (showst (tn,stn));
                 recur (tn,evhandle (stn,tol),tstep'))
        end
  in
   recur         
  end

val (tn,_) = driver (1000.0,make_stepper params,make_event params,hinterp)
                    (0.0,initial,tstep)

end
