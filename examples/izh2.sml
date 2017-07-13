(*
   Izhikevich Regular Spiking (RS) model with adaptive integration.
*)


structure IzhikevichRS2 =
struct

open RungeKutta

type IzhState = (real * real)

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

(* scaler :: real -> IzhState -> IzhState *)
fun scaler (s,(v,u)) = let open Real in (s*v, s*u) end

(* summer :: IzhState -> IzhState -> IzhState *)
fun summer ((v1,u1), (v2,u2)) = let open Real in (v1 + v2, u1 + u2) end


(* Parameters of the model *)

val params =
    {
      theta = 25.0,
      k1 = 0.04,
      k2 = 5.0,
      k3 = 140.0,
      a = 0.02,
      b = 0.2,
      c = ~65.0,
      d = 8.0,
      I = 10.0
    }


fun deriv {theta,k1,k2,k3,a,b,c,d,I}
	  (t,(v,u)) =
    let open Real
    in
        ( ((k1 * v * v) + (k2 * v) + k3 + (~ u) + I),
          (a * ((b * v) - u)))
    end

fun make_event (params as {theta,k1,k2,k3,a,b,c,d,I})  =
    let open Real 
    in
        (fn (v,u) =>
            ((*putStrLn ("make_event: v = " ^ (Real.toString v));*)
             (v - theta)),
         fn ((v,u),tstep) => 
            let
                val v' = c
                val u' = u + d
            in
                (v',u')
            end)
    end

datatype ('a, 'b) either = Left of 'a | Right of 'b

val tol = 1E~14
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


datatype 'a result = Next of 'a | Root of 'a


                                              
fun solver (stepper,evtest,hinterp) =
  let
      fun step (t,st,tstep) =
        let open Real
            val ((v,u),e as (ev,eu),ks) = stepper tstep (t,st)
        in
            case predictor tol (tstep,e) of
                Right tstep' => 
                (if (evtest (v,u) >= 0.0)
                 then (let
                          val tn      = t+tstep
                          val finterp = hinterp (tstep, ks, t, st)
                          val theta   = RootFind.brent tol (evtest o finterp) 0.0 1.0
                          val st'     = finterp theta
                      in
                          (putStrLn ("# event: theta = " ^ (showReal theta) ^ " tstep = " ^ (showReal tstep));
                           putStrLn ("# event: finterp theta = " ^ (showReal (#1(st'))));
                           Root (t+tstep*theta,st',tstep'))
                      end)
                 else Next (t+tstep,(v,u),tstep'))
              | Left tstep'  => 
                step (t,st,tstep')
        end
  in
      step
  end
    


val initial = (~65.0,~13.0)

val tstep = 0.5

val cerkdp: IzhState stepper3 = make_cerkdp()
val hinterp: IzhState hinterp = make_interp_cerkdp () (scaler,summer)
fun make_stepper (params) = cerkdp (scaler,summer,deriv params)


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
