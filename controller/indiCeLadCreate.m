function cel = indiCeLadCreate( aircraft, modes_cntrl_idx )
% Constant parameters for LAD control effectiveness

% Disclamer:
%   SPDX-License-Identifier: GPL-3.0-only
% 
%   Copyright (C) 2023 Yannic Beyer
%   Copyright (C) 2023 TU Braunschweig, Institute of Flight Guidance
% *************************************************************************


cntrl_input_idx = aircraft.wing_main.geometry.segments.control_input_index_local;

idx_min = min(cntrl_input_idx(2,cntrl_input_idx(2,:)>0));
idx_max = max(cntrl_input_idx(2,:));
num_lad_cntrl = idx_max-idx_min+1;
% control input mapping
cel.cntrlim = zeros(size(cntrl_input_idx,2),num_lad_cntrl);
for i = 1:size(cntrl_input_idx,2)
    if cntrl_input_idx(2,i) > 0
        cel.cntrlim(i,cntrl_input_idx(2,i)-idx_min+1) = 1;
    end
end

% local panel chord and area
cel.ci = aircraft.wing_main.geometry.ctrl_pt.c;
cel.si = aircraft.wing_main.geometry.ctrl_pt.c.*diff(aircraft.wing_main.geometry.line_25.pos(2,:));
cel.yi = aircraft.wing_main.geometry.ctrl_pt.pos(2,:);
cel.cos_sweep = cos(aircraft.wing_main.interim_results.sweep);

% normalized wing downwash
AIC = aircraft.wing_main.interim_results.AIC_b+aircraft.wing_main.interim_results.AIC_t;
spanwise_mapping = ...
    pinv( ...
        AIC ...
        * diag(aircraft.wing_main.state.geometry.ctrl_pt.c/aircraft.wing_main.params.b/2*(2*pi)) ...
        );
    
cel.spanwise_mapping = spanwise_mapping;
cel.m = aircraft.eom_flexible.body.m;

% derivative of structural mode acceleration w.r.t. local lift
cel.daccdl = ...
    aircraft.eom_flexible.structure_red.M_inv(6+modes_cntrl_idx,:) ...
    * aircraft.wing_main.aeroelasticity.T_sc(:,3:4:end) ...
    * spanwise_mapping;

cel.daccalldl = ...
    aircraft.eom_flexible.structure_red.M_inv(7:end,7:end) ...
    * aircraft.wing_main.aeroelasticity.T_sc(7:end,3:4:end) ...
    * spanwise_mapping;

% derivative of structural mode acceleration w.r.t. local pitching moment
cel.daccdm = ...
    aircraft.eom_flexible.structure_red.M_inv(6+modes_cntrl_idx,:) ...
    * aircraft.wing_main.aeroelasticity.T_sc(:,4:4:end) ...
    * -spanwise_mapping;

cel.M = aircraft.eom_flexible.structure_red.M;
cel.K = aircraft.eom_flexible.structure_red.K;
cel.T_sc = aircraft.wing_main.aeroelasticity.T_sc;
cel.spanwise_mapping = spanwise_mapping;

end
