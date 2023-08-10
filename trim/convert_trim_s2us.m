function [trimPoint] = convert_trim_s2us(trimPoint, sim_step_out, model_description, downwash)%, model_description) 
% NB: model_description corresponds to the UNSTEADY model
x_wing  = sim_step_out.x_unst_trim_main(:);
X_wing  = sim_step_out.X_unst_trim_main(:);
z_wing  = sim_step_out.z_unst_trim_main(:);
z2_wing = sim_step_out.z2_unst_trim_main(:);

x_htp   = sim_step_out.x_unst_trim_htp(:);
X_htp   = sim_step_out.X_unst_trim_htp(:);
z_htp   = sim_step_out.z_unst_trim_htp(:);
z2_htp  = sim_step_out.z2_unst_trim_htp(:);

x_vtp   = sim_step_out.x_unst_trim_vtp(:);
X_vtp   = sim_step_out.X_unst_trim_vtp(:);
z_vtp   = sim_step_out.z_unst_trim_vtp(:);
z2_vtp  = sim_step_out.z2_unst_trim_vtp(:);

R_Ab_fus        = sim_step_out.R_Ab_fus_trim(:);
alpha_ind_wing  = sim_step_out.alpha_ind_main_wing_trim(:);
alpha_ind_htp   = sim_step_out.alpha_ind_htp_trim(:);
alpha_ind_vtp   = sim_step_out.alpha_ind_vtp_trim(:);

trimPoint = Add_UnstAeroStates(trimPoint, x_wing, X_wing, z_wing, z2_wing, alpha_ind_wing,'wing');
trimPoint = Add_UnstAeroStates(trimPoint, x_htp, X_htp, z_htp, z2_htp, alpha_ind_htp,'htp');
trimPoint = Add_UnstAeroStates(trimPoint, x_vtp, X_vtp, z_vtp, z2_vtp, alpha_ind_vtp,'vtp');
trimPoint = Add_UnstAeroStates_Fuselage(trimPoint, R_Ab_fus);

trimPoint = init_downwash(trimPoint, downwash, sim_step_out.Gamma_mean); 

trimPoint = tpConvert(trimPoint, model_description);

end


function [alpha,beta, Vt] = Calc_AirData(states_init, model_description)
    states_list_full = [model_description.States.Continuous(:,1); model_description.States.Discrete(:,1)];    

    idx_q = dl2idx(states_list_full, {'qt_0';'qt_1';'qt_2';'qt_3'});
    q_bg = states_init(idx_q);
    EulerAngles  = q_bg2Phi( q_bg );
    
    idx_V = dl2idx(states_list_full, {'u_Kb';'v_Kb';'w_Kb'});
%     M_bg = Phi2DCM(EulerAngles);
    V_Ab = states_init(idx_V); %M_bg * V_Kg;
    [ alpha, beta ] = aeroAngles( V_Ab );
    
    Vt = norm(V_Ab,2);
    
end

function [trimPoint] = Add_UnstAeroStates(trimPoint, unst_aero_x,...
                        unst_aero_X, unst_aero_z, unst_aero_z2, alpha_ind, id )
    
    trimPoint = tpRmVal(trimPoint,'state',['alpha_ind_',id,'__']);
    
    trimPoint = tpAddVal(trimPoint,'continuous state',['x_aero_',id,'__'],unst_aero_x);
    trimPoint = tpAddVal(trimPoint,'continuous state',['X_',id,'__'],unst_aero_X);
    trimPoint = tpAddVal(trimPoint,'continuous state',['flap_aero_',id,'__'],unst_aero_z);
    trimPoint = tpAddVal(trimPoint,'continuous state',['act_aero_',id,'__'],unst_aero_z2);
    trimPoint = tpAddVal(trimPoint,'continuous state',['alpha_ind_',id,'__'],alpha_ind);
    trimPoint = tpAddVal(trimPoint,'continuous state',['tau_v_',id,'__'],zeros(size(alpha_ind)));
    
end

function [trimPoint] = Add_UnstAeroStates_Fuselage(trimPoint, R_Ab_fus)

    trimPoint = tpAddVal(trimPoint,'continuous state','R_Ab_fus__',R_Ab_fus);
    
end

