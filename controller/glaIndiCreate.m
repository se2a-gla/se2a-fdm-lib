function gla_indi = glaIndiCreate( aircraft, fp, varargin )
% GLA INDI controller paramters struct

% Disclamer:
%   SPDX-License-Identifier: GPL-3.0-only
% 
%   Copyright (C) 2022-2023 Yannic Beyer
%   Copyright (C) 2022 TU Braunschweig, Institute of Flight Guidance
% *************************************************************************

LAD_omega = aircraft.actuators.LAD.omega;
LAD_defl_max = aircraft.actuators.LAD.defl_max;
LAD_rate_max = aircraft.actuators.LAD.defl_rate_max;

% sensor filter
gla_indi.sflt.omega = 2.0 * LAD_omega;
gla_indi.sflt.d = 1;

% actuator boost factors
f_b_servo = 3;
f_b_aero = f_b_servo;
f_b_aero_max = 2*f_b_aero;

gla_indi.mcidx = [1,7];

for i = 1:length(varargin)
    if strcmp(varargin{i},'ModeControlIdx')
        gla_indi.mcidx = varargin{i+1};
    end
end

W_v_1 = 1e-9;
W_v_2 = ones(size(gla_indi.mcidx));

for i = 1:length(varargin)
    if strcmp(varargin{i},'SensOmega')
        gla_indi.sflt.omega(:) = varargin{i+1};
    elseif strcmp(varargin{i},'BoostServos')
        f_b_servo(:) = varargin{i+1};
    elseif strcmp(varargin{i},'BoostUnstAeroMin')
        f_b_aero(:) = varargin{i+1};
    elseif strcmp(varargin{i},'BoostUnstAeroMax')
        f_b_aero_max(:) = varargin{i+1};
    elseif strcmp(varargin{i},'WeightModes')
        W_v_2(:) = varargin{i+1};
    elseif strcmp(varargin{i},'WeightAcc')
        W_v_1(:) = varargin{i+1};
    end
end

f_b_aero(f_b_aero>f_b_aero_max) = f_b_aero_max;

num_cntrl_modes = length(gla_indi.mcidx);
% LAD control effectiveness
gla_indi.ce = indiCeLadCreate( aircraft, gla_indi.mcidx );

% unsteady flap time constant
atm     = isAtmosphere(fp.Altitude);
V       = fp.EAS * sqrt(1.225/atm.rho);
Ma      = V/atm.a;
V       = cos(aircraft.wing_main.interim_results.sweep)*V;
Ma      = cos(aircraft.wing_main.interim_results.sweep50)*Ma;
c       = aircraft.wing_main.geometry.ctrl_pt.c./cos(aircraft.wing_main.interim_results.sweep);
b_1     = 0.310;
b_2     = 0.312;
beta    = sqrtReal(1-powerFast(Ma,2));
beta2   = powerFast(beta,2);
Vc      = 2*V./c;
a_22    = -(b_1+b_2)*Vc.*beta2;

T_aero_vec  = -1 ./ (a_22/2);
T_aero = min(T_aero_vec);

cntrl_input_idx = aircraft.wing_main.geometry.segments.control_input_index_local;

idx_min = min(cntrl_input_idx(1,cntrl_input_idx(1,:)>0));
idx_max = max(cntrl_input_idx(1,:));
num_lad_cntrl = idx_max-idx_min+1;
% control input mapping
gla_indi.atc = zeros(1,num_lad_cntrl);
for i = idx_min:idx_max
    is_current_input = cntrl_input_idx(1,:) == i;
    gla_indi.atc(i) = T_aero_vec(is_current_input);
end

gla_indi.dtc = 2/LAD_omega;


% boost high-pass filter parameters
T1 = 1/LAD_omega;
bfa = f_b_aero * gla_indi.atc/T_aero;

gla_indi.boost_servo.b1 = f_b_servo;
gla_indi.boost_servo.b0 = f_b_servo/T1;
gla_indi.boost_servo.a1 = 1;
gla_indi.boost_servo.a0 = f_b_servo/T1;

gla_indi.boost_aero.b1 = bfa;
gla_indi.boost_aero.b0 = bfa./gla_indi.atc;
gla_indi.boost_aero.a1 = 1;
gla_indi.boost_aero.a0 = bfa./gla_indi.atc;

% desired closed-loop eigenfrequency
T_h = 1/f_b_servo*(2/LAD_omega) + 1/f_b_aero*T_aero + 2/gla_indi.sflt.omega;
mode_freq_1 = 1/3/T_h;
mode_freq = repmat(mode_freq_1,1,num_cntrl_modes); % -1 disables feedback control

% control allocation
gla_indi.ca = caLadParams( 1, -1, num_cntrl_modes+1, aircraft.wing_main.params.num_lads );
gla_indi.ca.W_v = diag([W_v_1,W_v_2]);
nl2 = aircraft.wing_main.params.num_lads/2;
W_u = [diff(aircraft.wing_main.geometry.line_25.pos(2,1:nl2+1)),diff(aircraft.wing_main.geometry.line_25.pos(2,end-nl2:end))];
W_u = W_u/mean(W_u);
gla_indi.ca.W_u = diag(W_u);
% gla_indi.ca.gamma = 100;
% gla_indi.ca.W_u = diag(ones(aircraft.wing_main.params.num_lads,1));
% gla_indi.ca.W_u = diag((W_u+2)/3);

% sample time
gla_indi.ts = 1/600;

% rate limit
T_delay = 1/f_b_servo*gla_indi.dtc + 1/f_b_aero*gla_indi.atc(1) + 2/gla_indi.sflt.omega;
Delta_t = T_delay;
u_rate_max = 1/LAD_defl_max * LAD_rate_max;
gla_indi.Delta_u_max = u_rate_max * Delta_t;

% deflection reference filter
gla_indi.ref.omega = 2/10;
gla_indi.ref.d = 1;

% dead band
gla_indi.udead = 0.0;

% LAD feedback gains
T_h = 1/f_b_servo*(2/LAD_omega) + 1/f_b_aero*T_aero + 2/gla_indi.sflt.omega;
gla_indi.k.pos = zeros(num_cntrl_modes,1);
gla_indi.k.vel = zeros(num_cntrl_modes,1);
gla_indi.k.acc = zeros(num_cntrl_modes,1);
for i = 1:num_cntrl_modes
    if mode_freq(i) == -1
        p = zeros(1,3);
    else
        p = repmat(-mode_freq(i),1,3);
    end
    k=ndiFeedbackGainPlace(p,T_h,false);
    gla_indi.k.pos(i) = k(1);
    gla_indi.k.vel(i) = k(2);
    gla_indi.k.acc(i) = k(3);
end

end
