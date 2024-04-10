% ** Parameters for wls control allocation (SE2A MR Bwd V4) **

% Disclamer:
%   SPDX-License-Identifier: GPL-2.0-only
% 
%   Copyright (C) 2022 Yannic Beyer
%   Copyright (C) 2022 TU Braunschweig, Institute of Flight Guidance
% *************************************************************************

% k is the number of control inputs
% m is the number of pseudo control inputs

% minimum control input kx1 vector
param.u_min = [-1;-1;-1];
% maximum control input kx1 vector
param.u_max = ones(3,1);
% desired control input kx1 vector
param.u_d = zeros(3,1);

% weighting mxm matrix of pseudo-control
param.W_v = diag([1,0.01]);
% weighting kxk matrix of the control input vector
param.W_u = diag([100,100,1]);
% weighting of pseudo-control vs. control input (scalar)
param.gamma = 10000;
% initial working set kx1 vector
param.W = zeros(3,1);
% maximum number of iterations (scalar)
param.i_max = 100;
