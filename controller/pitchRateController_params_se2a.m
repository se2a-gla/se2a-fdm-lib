% ** INDI attitude flight mode parameters (default) **

% Disclamer:
%   SPDX-License-Identifier: GPL-2.0-only
% 
%   Copyright (C) 2020-2022 Yannic Beyer
%   Copyright (C) 2022 TU Braunschweig, Institute of Flight Guidance
% *************************************************************************

param.sflt.omega = 90;
% damping ratio of the second order low pass filter
param.sflt.D = 1;

param.servo.D = 1;
param.servo.omega = 30;

% attitude controller
% param.atc = loadParams( 'cntrlAttiRedIndi_params_se2a' );

param.rm.T = 1;

k = ndiFeedbackGainPlace(-5*[1,1],2/(param.servo.omega)+2/(param.sflt.omega));
param.k.rate = k(1);
param.k.acc  = k(2);

% control effectiveness
param.cntrl_effect = loadParams( 'control_effectiveness_params_se2a' );
param.cntrl_effect.ny_du_red = param.cntrl_effect.ny_du_red([2,5],1:3);
param.cntrl_effect.ny_du_dt = param.cntrl_effect.ny_du_dt([2,5],1:3);

% control allocation
param.ca = loadParams( 'ca_pitch_params_se2a' );

% flight mode sample time, in s
param.sample_time = 1/600;

