  function xdot = dist_dyn(t,x);
%
% solve for the transient stage compositons in an ideal
% binary distillation column using ode45.
% 
% (c) 1997 B. Wayne Bequette - 24 Jan 1997
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
% [t,x] = ode45('dist_dyn',t0,tf,x0)
%
% where x0 is a vector of initial values for the liquid
% phase stage compositions (length(x0) = ns)
%
  global DIST_PAR PLOT_REFLUX TIME_REFLUX
%
% DIST_PAR is a vector of distillation column parameters
%          used by both dist_ss.m and dist_dyn.m
%
  if length(DIST_PAR) < 11;
    disp('not enough parameters given in DIST_PAR')
	disp(' ')
	disp('check to see that global DIST_PAR has been defined')
	return
  end
%
  alpha   = DIST_PAR(1);  % relative volatility (1.5)
  ns      = DIST_PAR(2);  % total number of stages (41)
  nf      = DIST_PAR(3);  % feed stage (21)
  feedi   = DIST_PAR(4);  % initial feed flowrate (1)
  zfeedi  = DIST_PAR(5);  % initial feed composition, light comp (0.5)
  qf      = DIST_PAR(6);  % feed quality (1 = sat'd liqd,
%                              0 = sat'd vapor) (1)
  refluxi = DIST_PAR(7);  % initial reflux flowrate (2.706)
  vapori  = DIST_PAR(8);  % initial reboiler vapor flowrate (3.206)
  md      = DIST_PAR(9);  % distillate molar hold-up (5)
  mb      = DIST_PAR(10); % bottoms molar hold-up (5)
  mt      = DIST_PAR(11); % stage molar hold-up (0.5)
%
  if length(DIST_PAR) >= 21;
   stepr   = DIST_PAR(12); % magnitude step in reflux (0)
   tstepr  = DIST_PAR(13); % time of reflux step change (0)
   stepv   = DIST_PAR(14); % magnitude step in vapor (0)
   tstepv  = DIST_PAR(15); % time of vapor step change (0)
   stepzf  = DIST_PAR(16); % magnitude of feed comp change (0)
   tstepzf = DIST_PAR(17); % time of feed comp change (0)
   stepf   = DIST_PAR(18); % magnitude of feed flow change (0)
   tstepf  = DIST_PAR(19); % time of feed flow change (0)
   step2r   = DIST_PAR(20); % magnitude of feed flow change (0)
   tstep2r  = DIST_PAR(21); % time of feed flow change (0)
  else
   stepr = 0; tstepr = 0; stepv = 0; tstepv = 0;
   stepzf = 0; tstepzf = 0; stepf = 0; tstepf = 0; step2r = 0; tstep2r = 0;
  end
%
% DIST_PAR(9:19) used by dist_dyn.m (distillation dynamics)
%  dist     = distillate product flowrate
%  lbot     = bottoms product flowrate
%  lr       = liquid flow in rectifying section (top)
%  ls       = liquid flow in stripping section (bottom)
%  vr       = vapor flow - rectifying sec (= vapor + feed*(1-qf)
%  vs       = vapor flow - stripping section (= vapor)
%  x(i)     = mole frac light component on stage i, liq
%  xdot(i)  = light component ith stage mat bal equation
%  y(i)     = mole frac light component on stage i, vap
%
%  check disturbances in reflux, vapor boil-up, feed composition
%        and feed flowrate
%
% IMPLEMENT TWO STEP CHANGES IN FORCING FUNCTION
% (reflux flowrate)
%
  if t < tstepr;
    reflux = refluxi;
  elseif (t >= tstepr) && (t < tstep2r)
    reflux = refluxi+stepr;
  elseif (t >= tstep2r)
    reflux = refluxi+step2r;
  end
 %
 % Need to do something like this if you want to plot the variable(s)
 % undergoing step changes.
 %
 TIME_REFLUX = [TIME_REFLUX' t]';
 PLOT_REFLUX = [PLOT_REFLUX' reflux]';
 % 
  if t < tstepv;
    vapor = vapori;
  else
    vapor = vapori + stepv;
  end
%
  if t < tstepzf;
    zfeed = zfeedi;
  else
    zfeed = zfeedi + stepzf;
  end
%
  if t < tstepf;
    feed = feedi;
  else
    feed = feedi + stepf;
  end
%
% rectifying and stripping section liquid flowrates
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
     xdot = zeros(ns,1);
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
      xdot(1)=(1/md)*(vr*y(2)-(dist+reflux)*x(1));
%
% rectifying (top) section
%
      for i=2:nf-1;
        xdot(i)=(1/mt)*(lr*x(i-1)+vr*y(i+1)-lr*x(i)-vr*y(i));
      end
%
% feed stage
%
      xdot(nf)  = (1/mt)*(lr*x(nf-1)+vs*y(nf+1)-ls*x(nf)-vr*y(nf)+feed*zfeed);
%
% stripping (bottom) section
%
      for i=nf+1:ns-1;
        xdot(i)=(1/mt)*(ls*x(i-1)+vs*y(i+1)-ls*x(i)-vs*y(i));
      end
%
% reboiler
%
      xdot(ns)=(1/mb)*(ls*x(ns-1)-lbot*x(ns)-vs*y(ns));
