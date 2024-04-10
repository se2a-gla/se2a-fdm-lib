function param = caLadParams(u_max,u_min,num_cntrl_modes,num_lads)
% ** Parameters for wls control allocation (SE2A MR Bwd V4) **
% Inputs:
%   k is the number of control inputs
%   m is the number of pseudo control inputs

% Disclamer:
%   SPDX-License-Identifier: GPL-2.0-only
% 
%   Copyright (C) 2022 Yannic Beyer
%   Copyright (C) 2022 TU Braunschweig, Institute of Flight Guidance
% *************************************************************************

% minimum control input kx1 vector
param.u_min = u_min*ones(num_lads,1);
% maximum control input kx1 vector
param.u_max = u_max*ones(num_lads,1);
% desired control input kx1 vector
param.u_d = zeros(num_lads,1);

% weighting mxm matrix of pseudo-control
param.W_v = diag([1e-9,ones(1,num_cntrl_modes)]);
% weighting kxk matrix of the control input vector
param.W_u = eye(num_lads);
% weighting of pseudo-control vs. control input (scalar)
param.gamma = 1e7;
% initial working set kx1 vector
param.W = zeros(num_lads,1);
% maximum number of iterations (scalar)
param.i_max = 100;
