  function f = dist_ss(x);
%
% solve for the steady-state stage compositons in an ideal
% binary distillation column using fsolve.
% 
% (c) 1993 B. Wayne Bequette - 21 june 93
% revised 31 Dec 96
%
% All flowrates are molar quantities.  Stages are numbered
% from the top down.  A total condenser is assumed.
% The overhead receiver is stage 1.  The partial reboiler
% is stage ns (the number of equilibrium "trays" is then
% ns-1). The column parameters should be specified in the
% DIST_PAR array.
% 
% to use this function, enter the following in the command
% window, or from a script file (after defining parameters
% in the DIST_PAR array:
%
% x = fsolve('dist_ss',x0)
%
% where x0 is a vector of initial guesses for the liquid
% phase stage compositions (length(x0) = ns)
%
  global DIST_PAR
%
% DIST_PAR is a vector of distillation column parameters
%          used by both dist_ss.m and dist_dyn.m
%
  if length(DIST_PAR) < 8;
    disp('not enough parameters given in DIST_PAR')
	disp(' ')
	disp('check to see that global DIST_PAR has been defined')
	return
  end
%
  alpha   = DIST_PAR(1);  % relative volatility (2.5)
  ns      = DIST_PAR(2);  % total number of stages (3)
  nf      = DIST_PAR(3);  % feed stage (2)
  feed    = DIST_PAR(4);  % feed flowrate (1)
  zfeed   = DIST_PAR(5);  % feed composition, light comp (0.5)
  qf      = DIST_PAR(6);  % feed quality (1 = sat'd liqd,
%                              0 = sat'd vapor) (1)
  reflux  = DIST_PAR(7);  % reflux flowrate (3)
  vapor   = DIST_PAR(8);  % reboiler vapor flowrate (3.5)
%
% DIST_PAR(9:19) used by dist_dyn.m (distillation dynamics)
%  dist     = distillate product flowrate
%  f(i)     = ith comp mat bal equation
%  lbot     = bottoms product flowrate
%  lr       = liquid flow in rectifying section (top)
%  ls       = liquid flow in stripping section (bottom)
%  vr       = vapor flow - rectifying sec (= vapor + feed*(1-qf)
%  vs       = vapor flow - stripping section (= vapor)
%  x(i)     = mole frac light component on stage i, liq
%  y(i)     = mole frac light component on stage i, vap
%
% rectifying and stripping section liquid flowrates
%
     lr    =  reflux;
     ls    =  reflux + feed*qf;
%
% rectifying and stripping section vapor flowrates
%
     vs    =  vapor;
     vr    =  vs +  feed*(1-qf);
%
% distillate and bottoms rates
%
     dist  = vr - reflux;
     lbot  =  ls - vs;
%
     if dist < 0
	   disp('error in specifications, distillate flow < 0')
	   return
	 end
	 if lbot < 0
	   disp('error in specifications, stripping section ')
	   disp(' ')
	   disp('liquid flowrate is negative')
	   return
	 end
%
% zero the function vector
%
     f = zeros(ns,1);
%
% calculate the equilibrium vapor compositions
%
      for i=1:ns;
        y(i)=(alpha*x(i))/(1.+(alpha-1.)*x(i));
      end
%
% material balances
%
% overhead receiver
%
      f(1)=(vr*y(2)-(dist+reflux)*x(1));
%
% rectifying (top) section
%
      for i=2:nf-1;
        f(i)=lr*x(i-1)+vr*y(i+1)-lr*x(i)-vr*y(i);
      end
%
% feed stage
%
      f(nf)  =  lr*x(nf-1)+vs*y(nf+1)-ls*x(nf)-vr*y(nf)+feed*zfeed;
%
% stripping (bottom) section
%
      for i=nf+1:ns-1;
        f(i)=ls*x(i-1)+vs*y(i+1)-ls*x(i)-vs*y(i);
      end
%
% reboiler
%
      f(ns)=(ls*x(ns-1)-lbot*x(ns)-vs*y(ns));
