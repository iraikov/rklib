
fun putStrLn str = 
    (TextIO.output (TextIO.stdOut, str);
     TextIO.output (TextIO.stdOut, "\n"))
    
fun putStr str = 
    (TextIO.output (TextIO.stdOut, str))
    
fun showReal n = 
    let open StringCvt
	open Real
    in
	(if n < 0.0 then "-" else "") ^ (fmt (FIX (SOME 12)) (abs n))
    end

fun printstate ({t,tstep,v,w,Istim,vk,vl,vca,gk,gl,gca,c,v1,v2,v3,v4,phi})=
    ( (showReal t) ^ " " ^ (showReal v)  ^ (showReal w) )
    
fun run (tmax,input as {t,tstep,v,w,Istim,vk,vl,vca,gk,gl,gca,c,v1,v2,v3,v4,phi}) =
  let val nstate = Model.morris_lecar input
  in putStrLn (printstate nstate); 
      if (#t nstate)  > tmax
     then (putStrLn "# All done!"; nstate)
     else (run (tmax,nstate))
  end

val initial = 
    {t     = 0.0,
     tstep = 0.001,
     v     = ~60.899,
     w     = 0.0149,
     Istim = 20.0,
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


val _ = run (100.0,initial)
