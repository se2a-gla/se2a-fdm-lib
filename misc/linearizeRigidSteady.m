function lin_sys = linearizeRigidSteady(model_name, model_description, trim_point, method)
% Create linear model of aircraft
% Warning: linmod method does not take into account states defined in
% lin_config.


if nargin<4 %Use jj_lin by default
    method = 'jj_lin';
end

%% Define list of states, inputs, and outputs
x_names = trim_point.mdl_names.States.Continuous; 
u_names = trim_point.mdl_names.MDL_Inputs;
y_names = trim_point.mdl_names.MDL_Outputs;

numLinStates    = size(x_names, 1);
numLinInputs    = size(u_names, 1);
numLinOutputs   = size(y_names,1);

%% Find corresponding indices of lists
numStates     = size(trim_point.State, 1);
numContStates = size(model_description.States.Continuous,1);

X_idxs = dl2idx(model_description.States.Continuous, x_names, 1);

U_idxs = dl2idx(model_description.MDL_Inputs, u_names, 1);

Xdot_idxs = X_idxs;

Y_idxs = dl2idx(model_description.MDL_Outputs, y_names, 1);

X_U_idxs    = [X_idxs; U_idxs + numStates];
Xdot_Y_idxs = [Xdot_idxs; Y_idxs+numContStates];

switch method
    case 'jj_lin'
    %% Compute jacobian

        jaco = jj_lin(model_name, [trim_point.State;trim_point.Input], numStates, X_U_idxs, Xdot_Y_idxs);

    %% Extract SS matrices

        A = jaco(1:numLinStates, 1:numLinStates);
        B = jaco(1:numLinStates, (numLinStates+1):(numLinStates+numLinInputs));
        C = jaco((numLinStates+1):(numLinStates+numLinOutputs), 1:numLinStates);
        D = jaco((numLinStates+1):(numLinStates+numLinOutputs), (numLinStates+1):(numLinStates+numLinInputs));
        
        lin_sys = ss(A,B,C,D);
        lin_sys.StateName    = x_names;
        lin_sys.OutputName   = y_names;
        lin_sys.InputName    = u_names;

    case {'linmod','linearize'}
        
        switch lin_config.method
            case 'linmod'
                lin_out = linmod(model_name, trim_point.State, trim_point.Input);
            case 'linearize'
                lin_out = linearize(model_name);
        end
        
        sys_tmp = ss(lin_out.a,lin_out.b, lin_out.c, lin_out.d,lin_out.Ts); %
        if lin_out.Ts > 0
            sys_tmp = d2c(sys_tmp);
        end
        
        states_mdl_cbd = [model_description.States.Continuous; model_description.States.Discrete];
        
        ii = 1;
        while ii <= length(lin_out.StateName)
            idx_mdl = find(strcmp(lin_out.StateName(ii),states_mdl_cbd(:,5)));
            if ~isempty(idx_mdl)
                n_states = length(idx_mdl);
                sys_tmp.StateName(ii:ii+n_states-1) = states_mdl_cbd(idx_mdl,1);
                ii = ii+n_states;
            else
                sys_tmp.StateName{ii} = ['xlin__', num2str(ii)];
                ii = ii+1;
            end
        end
        
        sys_tmp.OutputName = model_description.MDL_Outputs(:,1);
        sys_tmp.InputName = model_description.MDL_Inputs(:,1);    
        lin_sys = sys_tmp(y_names, u_names);
end


end

