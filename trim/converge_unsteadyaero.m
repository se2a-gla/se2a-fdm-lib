function [state_vec] = ...%[unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v, Delta_alpha, alpha_ind] = ...
                converge_unsteadyaero(states_init, inputs_init, model_description)
% Warning: development in progress, function is not generalized. 
states_list_full = [model_description.States.Continuous(:,1); model_description.States.Discrete(:,1)];

xyz_cg_c = evalin('base','config.xyz_cg_c');
atmosphere = evalin('base', 'envir.atmosphere');

downwash = evalin('base', 'downwash');

[alpha, beta, V] = Calc_AirData(states_init, model_description);
idx_omega = dl2idx(states_list_full, {'p_Kb';'q_Kb';'r_Kb'});
omega = states_init(idx_omega);

actuators_pos = zeros(14,1);
actuators_rate = zeros(14,1);

V_Kb_dt = zeros(3,1);
omega_dt = zeros(3,1);
structure_state = [zeros(6,1); states_init(44:67); zeros(6,1); states_init(68:91)];
structure_accel = zeros(30,1);

%% Wing Calculation

% Assign inputs
% wing_struct = wing_main;
% alpha = -0.0127;
% beta = 0;
% V = 150;
% omega = zeros(3,1);

wing_struct = evalin('base', 'wing_main');
wing_struct.config.is_circulation_iteration = 1;
n_panels = wing_struct.n_panel;

incidence = 0;
xyz_cg =  xyz_cg_c + wing_struct.geometry.origin;

V_ext_local = zeros(3,wing_struct.n_panel);
V_ext_local_dt = zeros(3, wing_struct.n_panel);

if ~wing_struct.config.is_unsteady
    unst_aero_x = zeros(8, wing_struct.n_panel);
    unst_aero_X = zeros(3, wing_struct.n_panel);
    unst_aero_z = zeros(2, wing_struct.n_panel);
    unst_aero_z2 = zeros(1,wing_struct.n_panel);
    tau_v = zeros(1,wing_struct.n_panel);
    Delta_alpha = zeros(1, wing_struct.n_panel);
    alpha_ind = zeros(1,wing_struct.n_panel);
else
    idx_x_aero = dl2idx(states_list_full, {'x_aero_wing__1'});
    idx_X = dl2idx(states_list_full, {'X_wing__1'});
    idx_z = dl2idx(states_list_full, {'flap_aero_wing__1'});
    idx_z2 = dl2idx(states_list_full, {'act_aero_wing__1'});
    idx_tau = dl2idx(states_list_full, {'tau_v_wing__1'});
%     idx_delta = dl2idx(states_list_full, {'Delta_alpha_wing__1'});
    idx_alphaind = dl2idx(states_list_full, {'alpha_ind_wing__1'});

    unst_aero_x =  reshape(states_init(idx_x_aero:idx_x_aero+n_panels*8-1), 8, n_panels);
    unst_aero_X =  reshape(states_init(idx_X:idx_X+n_panels*3-1), 3, n_panels);
    unst_aero_z =  reshape(states_init(idx_z:idx_z+n_panels*2-1), 2, n_panels);
    unst_aero_z2  =   reshape(states_init(idx_z2:idx_z2+n_panels-1), 1, n_panels);
    tau_v  =   reshape(states_init(idx_tau:idx_tau+n_panels-1), 1, n_panels);
%     Delta_alpha =   reshape(states_init(idx_delta:idx_delta+n_panels-1), 1, n_panels);
    alpha_ind  =   reshape(states_init(idx_alphaind:idx_alphaind+n_panels-1), 1, n_panels);
end

% incidence and induced velocity
wing_struct.state.geometry.ctrl_pt.local_incidence = wing_struct.geometry.ctrl_pt.local_incidence + incidence;

% compute wing state
n_iter_max = 3000;
dt = 0.00033;
eps_max = 1e-6;
% [unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v, Delta_alpha, alpha_ind, wing_data] =... 
%             Run_ConvergenceLoop(n_iter_max, dt, eps_max, ...
%             wing_struct, alpha, beta,...
%             V, omega, actuators_pos, actuators_rate, xyz_cg, V_Kb_dt, omega_dt,...
%             atmosphere, V_ext_local, V_ext_local_dt,...
%             structure_state, structure_accel, ...
%             unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v,...
%             Delta_alpha, alpha_ind);

[unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v, alpha_ind, wing_data] =... %Delta_alpha, alpha_ind, wing_data] =... 
            Run_UnsteadyTrimLoop(n_iter_max, dt, eps_max, ...
            wing_struct, alpha, beta,...
            V, omega, actuators_pos, actuators_rate, xyz_cg, V_Kb_dt, omega_dt,...
            atmosphere, V_ext_local, V_ext_local_dt,...
            structure_state, structure_accel, ...
            unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v,...
            alpha_ind);  
      
state_vec = Assign_UnstAeroStates_Wing(states_init, unst_aero_x, unst_aero_X,...
    unst_aero_z, unst_aero_z2, tau_v, alpha_ind, model_description, wing_struct.n_panel,  wing_struct.config.is_unsteady);

%% HTP Calculation
wing_struct = evalin('base', 'wing_htp');
wing_struct.config.is_circulation_iteration = 1;
n_panels = wing_struct.n_panel;

idx_htp = dl2idx(model_description.MDL_Inputs(:,1), {'de_htp'});
incidence = inputs_init(idx_htp);
xyz_cg =  xyz_cg_c + wing_struct.geometry.origin;

% Downwash calculation
V_ext_local = zeros(3,wing_struct.n_panel);
V_ext_local_dt = zeros(3, wing_struct.n_panel);
Gamma_mean = mean(wing_data.state.aero.circulation.Gamma(abs(wing_data.geometry.ctrl_pt.pos(2,:))<wing_struct.params.b));
assignin('base', 'Gamma_mean', Gamma_mean);
w_downwash = Gamma_mean * evalin('base', 'downwash.k_st');
V_ext_local(3,:) = w_downwash;

if ~wing_struct.config.is_unsteady
    unst_aero_x = zeros(8, wing_struct.n_panel);
    unst_aero_X = zeros(3, wing_struct.n_panel);
    unst_aero_z = zeros(2, wing_struct.n_panel);
    unst_aero_z2 = zeros(1,wing_struct.n_panel);
    tau_v = zeros(1,wing_struct.n_panel);
    
%     idx_delta    = dl2idx(states_list_full, {'Delta_alpha_htp__1'});
%     idx_alphaind = dl2idx(states_list_full, {'alpha_ind_htp__1'});
%     Delta_alpha  =   reshape(states_init(idx_delta:idx_delta+n_panels-1), 1, n_panels);
%     alpha_ind    =   reshape(states_init(idx_alphaind:idx_alphaind+n_panels-1), 1, n_panels);    
    Delta_alpha = zeros(1, wing_struct.n_panel);
    alpha_ind = zeros(1,wing_struct.n_panel);
else
    idx_x_aero = dl2idx(states_list_full, {'x_aero_htp__1'});
    idx_X = dl2idx(states_list_full, {'X_htp__1'});
    idx_z = dl2idx(states_list_full, {'flap_aero_htp__1'});
    idx_z2 = dl2idx(states_list_full, {'act_aero_htp__1'});
    idx_tau = dl2idx(states_list_full, {'tau_v_htp__1'});
%     idx_delta = dl2idx(states_list_full, {'Delta_alpha_htp__1'});
    idx_alphaind = dl2idx(states_list_full, {'alpha_ind_htp__1'});

    unst_aero_x  =  reshape(states_init(idx_x_aero:idx_x_aero+n_panels*8-1), 8, n_panels);
    unst_aero_X  =  reshape(states_init(idx_X:idx_X+n_panels*3-1), 3, n_panels);
    unst_aero_z  =  reshape(states_init(idx_z:idx_z+n_panels*2-1), 2, n_panels);
    unst_aero_z2 =   reshape(states_init(idx_z2:idx_z2+n_panels-1), 1, n_panels);
    tau_v        =   reshape(states_init(idx_tau:idx_tau+n_panels-1), 1, n_panels);
%     Delta_alpha  =   reshape(states_init(idx_delta:idx_delta+n_panels-1), 1, n_panels);
    alpha_ind    =   reshape(states_init(idx_alphaind:idx_alphaind+n_panels-1), 1, n_panels);
end

% incidence and induced velocity
wing_struct.state.geometry.ctrl_pt.local_incidence = wing_struct.geometry.ctrl_pt.local_incidence + incidence;

% compute wing state
n_iter_max = 1000;
dt = 0.001;
eps_max = 1e-6;
% [unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v, Delta_alpha, alpha_ind, ~] =... 
%             Run_ConvergenceLoop(n_iter_max, dt, eps_max, ...
%             wing_struct, alpha, beta,...
%             V, omega, actuators_pos, actuators_rate, xyz_cg, V_Kb_dt, omega_dt,...
%             atmosphere, V_ext_local, V_ext_local_dt,...
%             structure_state, structure_accel, ...
%             unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v,...
%             Delta_alpha, alpha_ind);
        
[unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v, alpha_ind, ~] =... 
            Run_UnsteadyTrimLoop(n_iter_max, dt, eps_max, ...
            wing_struct, alpha, beta,...
            V, omega, actuators_pos, actuators_rate, xyz_cg, V_Kb_dt, omega_dt,...
            atmosphere, V_ext_local, V_ext_local_dt,...
            structure_state, structure_accel, ...
            unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v,...
            alpha_ind);
        
state_vec = Assign_UnstAeroStates_HTP(state_vec, unst_aero_x, unst_aero_X,...
    unst_aero_z, unst_aero_z2, tau_v, alpha_ind, model_description, wing_struct.n_panel, wing_struct.config.is_unsteady);

end
%% get aerodynamic forces and moments

% XYZ_i_b         = wingGetLocalForce( wing_struct );
% LMN_i_b         = wingGetLocalMoment( wing_struct );
% M_profile_i_b   = wingGetLocalAirfoilMoment( wing_struct );
% XYZ_b           = wingGetGlobalForce( wing_struct );
% LMN_b           = wingGetGlobalMoment( wing_struct );
% num_iter        = wing_struct.state.aero.circulation.num_iter;
% wing_state      = wing_struct.state;

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

function [state_vec] = Assign_UnstAeroStates_Wing(states_init, unst_aero_x,...
                            unst_aero_X, unst_aero_z, unst_aero_z2, tau_v,...
                            alpha_ind, model_description, n_panels,  b_unst)
    state_vec = states_init;
    states_list_full = [model_description.States.Continuous(:,1); model_description.States.Discrete(:,1)];
    
    if b_unst
        idx_x_aero = dl2idx(states_list_full, {'x_aero_wing__1'});
        idx_X = dl2idx(states_list_full, {'X_wing__1'});
        idx_z = dl2idx(states_list_full, {'flap_aero_wing__1'});
        idx_z2 = dl2idx(states_list_full, {'act_aero_wing__1'});
        idx_tau = dl2idx(states_list_full, {'tau_v_wing__1'});
    
        state_vec(idx_x_aero:idx_x_aero+n_panels*8-1) = unst_aero_x;
        state_vec(idx_X:idx_X+n_panels*3-1) = unst_aero_X;
        state_vec(idx_z:idx_z+n_panels*2-1) = unst_aero_z;
        state_vec(idx_z2:idx_z2+n_panels-1) = unst_aero_z2;
        state_vec(idx_tau:idx_tau+n_panels-1) = tau_v; 
    end
    
%         idx_delta = dl2idx(states_list_full, {'Delta_alpha_wing__1'});
        idx_alphaind = dl2idx(states_list_full, {'alpha_ind_wing__1'});

%         state_vec(idx_delta:idx_delta+n_panels-1) = Delta_alpha;
        state_vec(idx_alphaind:idx_alphaind+n_panels-1) = alpha_ind;
%     end
end


function [state_vec] = Assign_UnstAeroStates_HTP(states_init, unst_aero_x,...
                            unst_aero_X, unst_aero_z, unst_aero_z2, tau_v,...
                            alpha_ind, model_description, n_panels,  b_unst)
    state_vec = states_init;
    states_list_full = [model_description.States.Continuous(:,1); model_description.States.Discrete(:,1)];
    
    if b_unst
        idx_x_aero = dl2idx(states_list_full, {'x_aero_htp__1'});
        idx_X = dl2idx(states_list_full, {'X_htp__1'});
        idx_z = dl2idx(states_list_full, {'flap_aero_htp__1'});
        idx_z2 = dl2idx(states_list_full, {'act_aero_htp__1'});
        idx_tau = dl2idx(states_list_full, {'tau_v_htp__1'});
        
        state_vec(idx_x_aero:idx_x_aero+n_panels*8-1) = unst_aero_x;
        state_vec(idx_X:idx_X+n_panels*3-1) = unst_aero_X;
        state_vec(idx_z:idx_z+n_panels*2-1) = unst_aero_z;
        state_vec(idx_z2:idx_z2+n_panels-1) = unst_aero_z2;
        state_vec(idx_tau:idx_tau+n_panels-1) = tau_v;
    end
    
%         idx_delta = dl2idx(states_list_full, {'Delta_alpha_htp__1'});
        idx_alphaind = dl2idx(states_list_full, {'alpha_ind_htp__1'});

%         state_vec(idx_delta:idx_delta+n_panels-1) = Delta_alpha;
        state_vec(idx_alphaind:idx_alphaind+n_panels-1) = alpha_ind;
%     end
end

function [unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v,  alpha_ind, wing_struct] =... 
            Run_ConvergenceLoop(n_iter_max, dt, eps_max, wing_struct, alpha, beta,...
            V, omega, actuators_pos, actuators_rate, xyz_cg, V_Kb_dt, omega_dt,...
            atmosphere, V_ext_local, V_ext_local_dt,...
            structure_state, structure_accel, ...
            unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v,...
             alpha_ind)
    iteration = 0;
    max_dt = 999;
    tic
    while (max_dt > eps_max) && (iteration < n_iter_max)
        iteration = iteration + 1;

        wing_struct = wingSetState( wing_struct, alpha, beta, V, omega, actuators_pos, actuators_rate, ...
            xyz_cg, 'accel', V_Kb_dt, 'angular_accel', omega_dt, 'atmosphere', ...
            atmosphere, 'wind', V_ext_local, V_ext_local_dt, 'structure_pos', ...
            structure_state(1:end/2), 'structure_vel', structure_state(end/2+1:end), ...
            'structure_accel', structure_accel, ...
            'unst_airfoil_state', unst_aero_x, 'dyn_stall_state', unst_aero_X, ...
            'unst_flap_state', unst_aero_z, 'unst_act2_state', unst_aero_z2, ...
            'tau_v', tau_v, 'alpha_ind', alpha_ind );

    %     unst_aero_x_dt = wing_struct.state.aero.unsteady.x_dt;
    %     unst_aero_X_dt = wing_struct.state.aero.unsteady.X_dt;
    %     unst_aero_z_dt = wing_struct.state.aero.unsteady.z_dt;
    %     unst_aero_z2_dt = wing_struct.state.aero.unsteady.z2_dt;
    %     tau_v_dt        = wing_struct.state.aero.unsteady.tau_v_dt;

        unst_aero_x_dt = wing_struct.state.aero.unsteady.x_dt;
        unst_aero_X_dt = wing_struct.state.aero.unsteady.X_dt;
        unst_aero_z_dt = wing_struct.state.aero.unsteady.z_dt;
        unst_aero_z2_dt = wing_struct.state.aero.unsteady.z2_dt;
        tau_v_dt        = wing_struct.state.aero.unsteady.tau_v_dt;

        unst_aero_x = unst_aero_x + dt*unst_aero_x_dt;
        unst_aero_X = unst_aero_X + dt*unst_aero_X_dt;
        unst_aero_z = unst_aero_z + dt*unst_aero_z_dt;
        unst_aero_z2 = unst_aero_z2 + dt*unst_aero_z2_dt;
        tau_v_dt     = tau_v + dt*tau_v_dt;

%         Delta_alpha_old = Delta_alpha;
        alpha_ind_old = alpha_ind;
%         Delta_alpha = wing_struct.state.aero.circulation.Delta_alpha;
        alpha_ind = wing_struct.state.aero.circulation.alpha_ind;
        
        max_dt_old = max_dt;
        
        [max_dt_ind, idx_max] = max(abs([unst_aero_x_dt; unst_aero_X_dt; unst_aero_z_dt;unst_aero_z2_dt; tau_v_dt]), [], 2);
        [max_dt] = max(abs([unst_aero_x_dt; unst_aero_X_dt; unst_aero_z_dt;unst_aero_z2_dt; tau_v_dt]), [], 'all');

    %     disp(['Iteration: ', num2str(iteration), ', max_dt: ', num2str(max_dt)])
    %     toc
    end
    disp(['Iteration: ', num2str(iteration), ', max_dt: ', num2str(max_dt)])
    toc
end

function [unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v,  alpha_ind, wing_struct] =... 
            Run_UnsteadyTrimLoop(n_iter_max, dt, eps_max, wing_struct, alpha, beta,...
            V, omega, actuators_pos, actuators_rate, xyz_cg, V_Kb_dt, omega_dt,...
            atmosphere, V_ext_local, V_ext_local_dt,...
            structure_state, structure_accel, ...
            unst_aero_x, unst_aero_X, unst_aero_z, unst_aero_z2, tau_v,...
            alpha_ind)
    iteration = 0;
    max_dt = 999;
    tic
    while iteration< 3 %(max_dt > eps_max) && (iteration < n_iter_max)
        iteration = iteration + 1;

        wing_struct = wingSetState( wing_struct, alpha, beta, V, omega, actuators_pos, actuators_rate, ...
            xyz_cg, 'accel', V_Kb_dt, 'angular_accel', omega_dt, 'atmosphere', ...
            atmosphere, 'wind', V_ext_local, V_ext_local_dt, 'structure_pos', ...
            structure_state(1:end/2), 'structure_vel', structure_state(end/2+1:end), ...
            'structure_accel', structure_accel, ...
            'unst_airfoil_state', unst_aero_x, 'dyn_stall_state', unst_aero_X, ...
            'unst_flap_state', unst_aero_z, 'unst_act2_state', unst_aero_z2, ...
            'tau_v', tau_v, 'alpha_ind', alpha_ind );

        [x, X, z, z2] = wingStateGetUnstAeroTrimOp( wing_struct.state, wing_struct.airfoil, wing_struct.config );

        unst_aero_x  = x;
        unst_aero_X  = X;
        unst_aero_z  = z;
        unst_aero_z2 = z2';    

        wing_struct.state.aero.unsteady.x = x;
        wing_struct.state.aero.unsteady.X = X;
        wing_struct.state.aero.unsteady.z = z;
        wing_struct.state.aero.unsteady.z2 = z2;    
        
    %     unst_aero_x_dt = wing_struct.state.aero.unsteady.x_dt;
    %     unst_aero_X_dt = wing_struct.state.aero.unsteady.X_dt;
    %     unst_aero_z_dt = wing_struct.state.aero.unsteady.z_dt;
    %     unst_aero_z2_dt = wing_struct.state.aero.unsteady.z2_dt;
    %     tau_v_dt        = wing_struct.state.aero.unsteady.tau_v_dt;

        unst_aero_x_dt = wing_struct.state.aero.unsteady.x_dt;
        unst_aero_X_dt = wing_struct.state.aero.unsteady.X_dt;
        unst_aero_z_dt = wing_struct.state.aero.unsteady.z_dt;
        unst_aero_z2_dt = wing_struct.state.aero.unsteady.z2_dt;
        tau_v_dt        = wing_struct.state.aero.unsteady.tau_v_dt;

%         unst_aero_x = unst_aero_x + dt*unst_aero_x_dt;
%         unst_aero_X = unst_aero_X + dt*unst_aero_X_dt;
%         unst_aero_z = unst_aero_z + dt*unst_aero_z_dt;
%         unst_aero_z2 = unst_aero_z2 + dt*unst_aero_z2_dt;
%         tau_v_dt     = tau_v + dt*tau_v_dt;

%         Delta_alpha_old = Delta_alpha;
        alpha_ind_old = alpha_ind;
%         Delta_alpha = wing_struct.state.aero.circulation.Delta_alpha;
        alpha_ind = wing_struct.state.aero.circulation.alpha_ind;
        
        tau_v = wing_struct.state.aero.unsteady.tau_v_dt;
%         
%         max_dt_old = max_dt;
%         
%         [max_dt_ind, idx_max] = max(abs([unst_aero_x_dt; unst_aero_X_dt; unst_aero_z_dt;unst_aero_z2_dt; tau_v_dt]), [], 2);
        [max_dt] = max(abs([unst_aero_x_dt; unst_aero_X_dt; unst_aero_z_dt;unst_aero_z2_dt; tau_v_dt]), [], 'all');
%         del_Da = mean(abs(Delta_alpha-Delta_alpha_old));
        del_ai = mean(abs(alpha_ind - alpha_ind_old));
    %     disp(['Iteration: ', num2str(iteration), ', max_dt: ', num2str(max_dt)])
    %     toc
    end
    disp(['Iteration: ', num2str(iteration), ', max_dt: ', num2str(max_dt),...
        ', alpha_ind: ', num2str(del_ai)])
    toc
end