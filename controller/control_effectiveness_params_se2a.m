% ** Parameters for control effectiveness (default) **

% Disclamer:
%   SPDX-License-Identifier: GPL-2.0-only
% 
%   Copyright (C) 2020-2022 Yannic Beyer
%   Copyright (C) 2022 TU Braunschweig, Institute of Flight Guidance
% *************************************************************************

% k is the number of control inputs
% m is the number of pseudo control inputs

% derivative of the pseudo control kx1 vector ny w.r.t. the control input
% mx1 vector u (usually actuator dynamics are reduced) at the trim
% point (jacobi kxm matrix)
param.ny_du_red = [ ...
    1.4296   -1.4296    0.0000   -1.0864         0; ...
   -0.0739   -0.0739   -1.1509   -0.0531         0; ...
    0.0321   -0.0321    0.0000    0.7502         0; ...
    0.0520    0.0520    0.2051         0    2.4938; ...
   -1.7313   -1.7313   -3.9215         0    2.4938 ...
    ];

% derivative of the pseudo control kx1 vector ny w.r.t. the time derivative
% of the control input mx1 vector u_dt at the trim point (jacobi kxm 
% matrix)
param.ny_du_dt = zeros(5,5);

% trim equivalent airspeed
param.V_trim = 100;
