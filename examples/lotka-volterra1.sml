
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

structure LotkaVolterra1 =
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

val rk4b: LVState stepper1 = make_rk4b()
val stepper = rk4b (scaler,summer,deriv)

val step = 0.01

fun putStrLn str =
    (TextIO.output (TextIO.stdOut, str);
     TextIO.output (TextIO.stdOut, "\n"))

fun putStr str =
    (TextIO.output (TextIO.stdOut, str))

fun showReal n = Real.toString n

fun showst (t, (y1, y2)) = String.concat [(showReal t), " ", (showReal y1) , " ", (showReal y2) ]

fun solver (t,st) =
  let val (tn,stn) = stepper step (t,st)
  in putStrLn (showst (t,st));
     if t > 40.00001
     then putStrLn "# All done!"
     else solver (tn,stn)
  end


val _ = solver (0.0,(2.0,0.3))

end

