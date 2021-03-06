(* Module for testing Runge-Kutta integration of ODEs
   Based on code by  Uwe Hollerbach <uh@alumni.caltech.edu>
*)

structure RungeKuttaTest =
struct

open RungeKutta

val summer = Real.+
val scaler = Real.*

infix 7 */
infix 6 +/
infix //



val con = ~0.4
fun deriv (t,y) = con*y
val t0 = 0.0
val y0 = 1.75
fun exact t = y0*Real.Math.exp(con*(t - t0))

fun putStr str =
    (TextIO.output (TextIO.stdOut, str))

fun putStrLn str =
    (TextIO.output (TextIO.stdOut, str);
     TextIO.output (TextIO.stdOut, "\n"))

fun showReal n = Real.toString n

fun showst (t, y) = String.concat [(showReal y), "\t", (showReal (y - (exact t)))]

fun gen_soln1 (integrator,h,t,st) =
  let 
      val stn = integrator (t,st)
      val tn  = Real.+(t,h)
  in 
      if t >= 5.0
      then putStrLn (showst (tn,stn))
      else gen_soln1 (integrator,h,tn,stn)
  end

fun gen_soln2 (integrator,h,t,st) =
  let 
      val (stn,en) = integrator (t,st)
      val tn       = Real.+(t,h)
  in 
      if t >= 5.0
      then putStrLn (showst (tn,stn))
      else gen_soln2 (integrator,h,tn,stn)
  end

fun gen_soln2_fsal (integrator,h,t,st,k) =
  let 
      val (stn,en,kn) = integrator (t,st,k)
      val tn       = Real.+(t,h)
  in 
      if t >= 5.0
      then putStrLn (showst (tn,stn))
      else gen_soln2_fsal (integrator,h,tn,stn,SOME kn)
  end

fun do_case1 integrator n =
  let 
      val h = if n < 0 then Real.Math.pow (2.0,Real.fromInt (~n)) 
	      else Real.Math.pow (0.5,Real.fromInt n)
      val sep = if n <= 4 then "\t\t" else "\t"
  in
      putStr (String.concat [(showReal h), sep]);
      gen_soln1 (integrator h,h,t0,y0)
  end

fun solver1 (integrator,stats) =
  (putStrLn stats;
   putStrLn "# step yf delta";
   List.app (do_case1 (integrator (scaler,summer,deriv)))
	    (List.tabulate (18, fn x => x - 2));
   putStrLn "# All done!\n")

fun do_case2 integrator n =
  let 
      val h = if n < 0 then Real.Math.pow (2.0,Real.fromInt (~n)) 
	      else Real.Math.pow (0.5,Real.fromInt n)
      val sep = if n <= 4 then "\t\t" else "\t"
  in
      putStr (String.concat [(showReal h), sep]);
      gen_soln2 (integrator h,h,t0,y0)
  end

fun do_case2_fsal integrator n =
  let 
      val h = if n < 0 then Real.Math.pow (2.0,Real.fromInt (~n)) 
	      else Real.Math.pow (0.5,Real.fromInt n)
      val sep = if n <= 4 then "\t\t" else "\t"
  in
      putStr (String.concat [(showReal h), sep]);
      gen_soln2_fsal (integrator h,h,t0,y0,NONE)
  end

fun solver2 (integrator,stats) =
  (putStrLn stats;
   putStrLn "# step yf delta";
   List.app (do_case2 (integrator (scaler,summer,deriv)))
	    (List.tabulate (18, fn x => x - 2));
   putStrLn "# All done!\n")

fun solver2_fsal (integrator,stats) =
  (putStrLn stats;
   putStrLn "# step yf delta";
   List.app (do_case2_fsal (integrator (scaler,summer,deriv)))
	    (List.tabulate (19, fn x => x - 2));
   putStrLn "# All done!\n")

fun gen_soln3 (integrator,interp,h,t,st) =
  let 
      val (stn,en,ks) = integrator (t,st)
      val tn       = Real.+(t,h)
  in 
      if t >= 5.0
      then (putStr (showst (tn,stn));
            putStr ("\t" ^ (showReal en));
            putStrLn ("\t" ^ (showReal (interp (h,ks,t,st) 1.0))))
      else gen_soln3 (integrator,interp,h,tn,stn)
  end

fun do_case3 integrator interp n =
  let 
      val h = if n < 0 then Real.Math.pow (2.0,Real.fromInt (~n)) 
	      else Real.Math.pow (0.5,Real.fromInt n)
      val sep = if n <= 4 then "\t\t" else "\t"
  in
      putStr (String.concat [(showReal h), sep]);
      gen_soln3 (integrator h,interp,h,t0,y0)
  end

fun solver3 (integrator,interp,stats) =
  (putStrLn stats;
   putStrLn "# step yf delta err uf";
   List.app (do_case3 (integrator (scaler,summer,deriv)) (interp (scaler,summer)))
	    (List.tabulate (18, fn x => x - 2));
   putStrLn "# All done!\n")


val rkfe: real stepper1 = make_rkfe()
val rk3:  real stepper1 = make_rk3()
val rk4a: real stepper1 = make_rk4a()
val rk4b: real stepper1 = make_rk4b()

val rkhe:  real stepper2 = make_rkhe()
val rkbs:  real stepper2_fsal = make_rkbs()
val rkoz3:  real stepper2_fsal = make_rkoz3()
val rkn34: real stepper2 = make_rkn34()
val rkf45: real stepper2 = make_rkf45()
val rkck:  real stepper2 = make_rkck()
val rkoz4:  real stepper2_fsal = make_rkoz4()
val rkdp:  real stepper2_fsal = make_rkdp()
val rkdpb: real stepper2_fsal = make_rkdpb()
val rkf78: real stepper2 = make_rkf78()
val rkv65: real stepper2 = make_rkv65()

val rkhe_aux:  real stepper1  = make_rkhe_aux()
val rkbs_aux:  real stepper1  = make_rkbs_aux()
val rkoz3_aux: real stepper1  = make_rkoz3_aux()
val rkn34_aux: real stepper1  = make_rkn34_aux()
val rkf45_aux: real stepper1  = make_rkf45_aux()
val rkck_aux:  real stepper1  = make_rkck_aux()
val rkoz4_aux: real stepper1  = make_rkoz4_aux()
val rkdp_aux:  real stepper1  = make_rkdp_aux()
val rkdpb_aux:  real stepper1 = make_rkdpb_aux()
val rkf78_aux: real stepper1  = make_rkf78_aux()
val rkv65_aux: real stepper1  = make_rkv65_aux()

val cerkdp:  real stepper3 = make_cerkdp()
val cerkoz3:  real stepper3 = make_cerkoz3()
val cerkoz4:  real stepper3 = make_cerkoz4()

val interp_cerkdp  = make_interp_cerkdp()
val interp_cerkoz3 = make_interp_cerkoz3()
val interp_cerkoz4 = make_interp_cerkoz4()


                                          
fun run() =
 (putStrLn "#### Non-Adaptive Solvers";
  List.app solver1 [(rkfe, show_rkfe),
                    (rk3,  show_rk3),
                    (rk4a, show_rk4a),
                    (rk4b, show_rk4b)];
  putStrLn "#### Adaptive Solvers";
  List.app solver2 [(rkhe, show_rkhe),
                    (rkn34, show_rkn34),
                    (rkf45, show_rkf45),
                    (rkck, show_rkck),
                    (rkf78, show_rkf78),
                    (rkv65, show_rkv65)];
  putStrLn "#### Adaptive Solvers (FSAL)";
  List.app solver2_fsal [(rkbs, show_rkbs),
                         (rkoz3, show_rkoz3),
                         (rkoz4, show_rkoz4),
                         (rkdp, show_rkdp),
                         (rkdpb, show_rkdpb)];
  putStrLn "#### Continuous Solvers";
  List.app solver3 [(cerkoz3, interp_cerkoz3, show_cerkoz3)];
  List.app solver3 [(cerkoz4, interp_cerkoz4, show_cerkoz4)];
  List.app solver3 [(cerkdp, interp_cerkdp, show_cerkdp)];
  putStrLn "#### Auxiliary Solvers: Error Estimators from Adaptives";
  List.app solver1 [(rkhe_aux, show_rkhe_aux),
		    (rkbs_aux, show_rkbs_aux),
		    (rkoz3_aux, show_rkoz3_aux),
		    (rkn34_aux, show_rkn34_aux),
		    (rkf45_aux, show_rkf45_aux),
		    (rkck_aux, show_rkck_aux),
		    (rkoz4_aux, show_rkoz4_aux),
		    (rkdp_aux, show_rkdp_aux),
		    (rkdpb_aux, show_rkdpb_aux),
		    (rkf78_aux, show_rkf78_aux),
                    (rkv65_aux, show_rkv65_aux)]
  )

val _ = run()
end
