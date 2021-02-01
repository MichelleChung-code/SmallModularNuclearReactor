function dxdt = pid_odefun(time,x)
global TIME_PLOT Q_PLOT FIN_PLOT FOUT_PLOT TIN_PLOT VSP_PLOT TSP_PLOT
global FSS QSS VSS TSS TINSS RHO CP TREF KC1 KC2 TI2
global TIN_STEPTIME TIN_STEP TSP_STEPTIME1 TSP_STEP1 TSP_STEPTIME2 TSP_STEP2
global FIN_STEPTIME FIN_STEP Q_MIN Q_MAX FOUT_MIN FOUT_MAX
%
v = x(1); t = x(2); integ = x(3);
%
fbias = FSS; qbias = QSS; vsp = VSS; 
%
% Disturbances
%
if time < TIN_STEPTIME
    tin = TINSS;
else
    tin = TINSS+TIN_STEP;
end
%
if time < FIN_STEPTIME
    fin = FSS;
else
    fin = FSS+FIN_STEP;
end
%
if time < TSP_STEPTIME1
    tsp = TSS;
elseif (time >= TSP_STEPTIME1) && (time < TSP_STEPTIME2)
    tsp = TSS+TSP_STEP1;
else
    tsp = TSS+TSP_STEP2;
end
%
error1 = vsp-v; error2 = tsp-t;
%
% P-only level control (direct-acting, KC1 < 0)
%
fout = fbias+KC1*error1;
%
% PI temperature control (reverse-acting, KC2 > 0)
%
q = qbias+KC2*(error2+1/TI2*integ);
%
% Clamp manipulated variables
%
fout = min(fout,FOUT_MAX); fout = max(fout,FOUT_MIN);
q = min(q,Q_MAX); q = max(q,Q_MIN);
%
dxdt(1) = fin-fout;
dxdt(2) = fin/v*(tin-TREF)-fout/v*(t-TREF)+1/(RHO*v*CP)*q;
dxdt(3) = error2;
%
TIME_PLOT = [TIME_PLOT' time]'; Q_PLOT = [Q_PLOT' q]'; 
FIN_PLOT = [FIN_PLOT' fin]'; FOUT_PLOT = [FOUT_PLOT' fout]'; 
TIN_PLOT = [TIN_PLOT' tin]'; VSP_PLOT = [VSP_PLOT' vsp]'; 
TSP_PLOT = [TSP_PLOT' tsp]';  
%
dxdt = dxdt';



