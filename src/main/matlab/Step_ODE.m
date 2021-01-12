
% Function created by Periasamy Vijay
% Date 27/11/2017
% ------------------------------
% Description of input arguments:
% ------------------------------
% fhan - Function handle for the differential equation function
% Solver - String name of the ODE solver
% t_s - Step time
% t_t - Total simulation time
% Val_ini - Vector of initial values (before step) for the inputs or parameters
% Val_ini - Vector of final values (after step) for the inputs or parameters
% ini - Intial value vector for the differential equations
% ------------------------------
% Taken from: https://www.mathworks.com/matlabcentral/fileexchange/65201-step-response-of-a-non-linear-differential-equation-system
%
% Copyright (c) 2017, Vijay Periasamy
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


function [t,y] = Step_ODE(fhan, Solver, t_s, t_t, Val_ini, Val_fin, ini)

% before the step is applied

y0 = ini;

switch Solver
    case 'ode15s'
        
        tspan = [0 t_s];
        [t1,y1] = ode15s(@(t,y) fhan(t,y,Val_ini), tspan, y0);
        trows = size(t1,1); % number of rows of the time vector

% after the step is applied

        tspan = [0 t_t-t_s];
        y0 = y1(trows,:);
        [t2,y2] = ode15s(@(t,y) fhan(t,y,Val_fin), tspan, y0);
     
    case 'ode45'
        
        tspan = [0 t_s];
        [t1,y1] = ode45(@(t,y) fhan(t,y,Val_ini), tspan, y0);
        trows = size(t1,1); % number of rows of the time vector

% after the step is applied

        tspan = [0 t_t-t_s];
        y0 = y1(trows,:);
        [t2,y2] = ode45(@(t,y) fhan(t,y,Val_fin), tspan, y0);
        
    case 'ode23'
        
        tspan = [0 t_s];
        [t1,y1] = ode23(@(t,y) fhan(t,y,Val_ini), tspan, y0);
        trows = size(t1,1); % number of rows of the time vector

% after the step is applied

        tspan = [0 t_t-t_s];
        y0 = y1(trows,:);
        [t2,y2] = ode23(@(t,y) fhan(t,y,Val_fin), tspan, y0);
        
    case 'ode113'
        
        tspan = [0 t_s];
        [t1,y1] = ode113(@(t,y) fhan(t,y,Val_ini), tspan, y0);
        trows = size(t1,1); % number of rows of the time vector

% after the step is applied

        tspan = [0 t_t-t_s];
        y0 = y1(trows,:);
        [t2,y2] = ode113(@(t,y) fhan(t,y,Val_fin), tspan, y0);
        
    case 'ode23s'
        
        tspan = [0 t_s];
        [t1,y1] = ode23s(@(t,y) fhan(t,y,Val_ini), tspan, y0);
        trows = size(t1,1); % number of rows of the time vector

% after the step is applied

        tspan = [0 t_t-t_s];
        y0 = y1(trows,:);
        [t2,y2] = ode23s(@(t,y) fhan(t,y,Val_fin), tspan, y0);
        
    case 'ode23t'
        
        tspan = [0 t_s];
        [t1,y1] = ode23t(@(t,y) fhan(t,y,Val_ini), tspan, y0);
        trows = size(t1,1); % number of rows of the time vector

% after the step is applied

        tspan = [0 t_t-t_s];
        y0 = y1(trows,:);
        [t2,y2] = ode23t(@(t,y) fhan(t,y,Val_fin), tspan, y0);
        
    case 'ode23tb'
        tstep = 0.1; % todo pass this in so that it can apply to the other cases too
        tspan = 0:tstep:t_s;
        [t1,y1] = ode23tb(@(t,y) fhan(t,y,Val_ini), tspan, y0);
        trows = size(t1,1); % number of rows of the time vector

% after the step is applied

        tspan = 0:tstep:t_t-t_s;
        y0 = y1(trows,:);
        [t2,y2] = ode23tb(@(t,y) fhan(t,y,Val_fin), tspan, y0);
        
    otherwise
        
        disp('Error in the Solver name input argument')
end


t2 = t2 + t_s;  % time offseting

% concatenation to obtain final output

t = vertcat(t1,t2);
y = vertcat(y1,y2);


end

