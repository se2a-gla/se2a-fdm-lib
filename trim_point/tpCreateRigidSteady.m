function [trimPoint] = tpCreateRigidSteady(modelName, model_description, tp_task)

Y_Reqts     = tp_task.Y_Reqts;
Xdot_Reqts  = tp_task.Xdot_Reqts;
XCont_Var   = tp_task.XCont_Var;
XDisc_Var   = tp_task.XDisc_Var;
U_Var       = tp_task.U_Var;

%% Obtain initial state and model sizes:

[modelSizes, x_0] = feval(modelName,  [], [], [], 'sizes'); % See documentation for 'model' function.

numContStates   = modelSizes(1);
numDiscStates   = modelSizes(2);
numStates       = modelSizes(1) + modelSizes(2);
numOutputs      = modelSizes(3);
numInputs       = modelSizes(4);

%% Process the trim specifications
Y_Reqts    = Process_TrimSpec(Y_Reqts, model_description.MDL_Outputs);
Xdot_Reqts = Process_TrimSpec(Xdot_Reqts, model_description.States.Continuous);

XCont_Var = Process_TrimSpec(XCont_Var, model_description.States.Continuous);
XDisc_Var = Process_TrimSpec(XDisc_Var, model_description.States.Discrete);
U_Var     = Process_TrimSpec(U_Var, model_description.MDL_Inputs);

%% Check number of trim requirements vs variables
numReqts = size(Y_Reqts,1) + size(Xdot_Reqts,1);
numVar   = size(XCont_Var,1) + size(XDisc_Var,1) + size(U_Var, 1);

if (numVar - numReqts) ~= 0
    error('Number of trim variables/trim requirements do not match!');
end

%% Build up full vectors

%Find indices in model description lists
Y_Reqts_Idxs    = dl2idx(model_description.MDL_Outputs, Y_Reqts(:,1));

Xdot_Reqts_Idxs = dl2idx(model_description.States.Continuous, Xdot_Reqts(:,1));

XCont_Var_Idxs = dl2idx(model_description.States.Continuous,XCont_Var(:,1));
% XDisc_Var_Idxs = dl2idx(model_description.States.Discrete, XDisc_Var(:,1));
XDisc_Var_Idxs = []; % Nothing has been specified for discrete states, so this is left empty.
X_Var_Idxs     = [XCont_Var_Idxs; (XDisc_Var_Idxs + numContStates)]; % Combined continuous and discrete

U_Var_Idxs     = dl2idx(model_description.MDL_Inputs, U_Var(:,1));


Y_Names     = model_description.MDL_Outputs(:,1);
Xdot_Names  = strcat(model_description.States.Continuous(:,1), '_dot');
X_Names     = [model_description.States.Continuous(:,1); model_description.States.Discrete(:,1)];
U_Names     = model_description.MDL_Inputs(:,1);

%   Values of outputs/state derivatives which are not specified in the 
%   requirements are not used, so set to zero by default
y_0 = zeros(numOutputs,1);
xdot_0 = zeros(numContStates,1);

u_0 = zeros(numInputs,1); % Inputs are assumed set to 0 unless used for trimming.

y_0(Y_Reqts_Idxs)       = cell2mat(Y_Reqts(:,2));
xdot_0(Xdot_Reqts_Idxs) = cell2mat(Xdot_Reqts(:,2)); 
u_0(U_Var_Idxs)         = cell2mat(U_Var(:,2));
x_0(XCont_Var_Idxs)     = cell2mat(XCont_Var(:,2));
% x_0(XDisc_Var_Idxs + numContStates) = cell2mat(XDisc_Var(:,2));

% x_0 = converge_unsteadyaero(x_0, u_0, model_description);

%% Run trim
trim_opts = [];
trim_opts.CompileFlag = 1;
trim_opts.n_bisec_max = 1000;
trim_opts.n_iter_max = 1000;
trim_opts.cost_tbg = 2e-13;

[x_tr, u_tr, d_tr, y_tr, trim_fail, trim_info] = jj_trim(modelName, x_0, u_0, xdot_0, y_0, X_Var_Idxs, U_Var_Idxs, Xdot_Reqts_Idxs, Y_Reqts_Idxs, ...
        X_Names, U_Names, Xdot_Names, Y_Names, [],[],1e-6*(1 + abs(x_0)),1e-6*(1 + abs(u_0)),trim_opts);
    
if trim_fail
    warning('Trim unsuccessful!');
end


% assignin('base', 'x_tr_orig', x_tr);
% x_tr_rc = converge_unsteadyaero(x_tr, u_tr, model_description);
% delta_x = sum(abs(x_tr - x_tr_rc))
% x_tr = x_tr_rc; 


%% Set up trim point

trimPoint.State = x_tr;
trimPoint.Input = u_tr;
trimPoint.Output = y_tr;
trimPoint.Deriv = d_tr;

trimPoint.mdl_names = mdlGetMdlNames(model_description);

trim_info = rmfield(trim_info, {'jaco', 'cost', 'StepType', 'BisecCounter', 'LimitedStep', 'x', 'u','d','y'});
trimPoint.info = trim_info;
trimPoint.Task = tp_task;

end

function trimSpec_out = Process_TrimSpec(trimSpec_in, model_desc_sub)
       
    if isempty(trimSpec_in)
        numberedVarIdx = [];
    else   
        numberedVarIdx = find(strncmp('_', trimSpec_in(:,1),1)); % Find 
    end
    
    trimSpec_out = trimSpec_in;
    trimSpec_tmp = []; % To contain rows to be added to trimSpec
    
    rowsToBeDeleted = []; %Array of rows to be removed from trimSpec_out later on. 
    
    for in = 1:size(numberedVarIdx,1) % Loop through the specified 'numbered variables' which have been found
        nVarName  = trimSpec_in{numberedVarIdx(in),1}; % 
        nVarValue = trimSpec_in{numberedVarIdx(in),2};
        
        nVarMDLIdxs = dl2idx_single(model_desc_sub, [nVarName(2:end), '__'],1);
        
        rowsToBeDeleted = [rowsToBeDeleted; numberedVarIdx(in)];
        
        for jn = 1:size(nVarMDLIdxs,1) % Loop through the set of variables connected to the 'numbered variable'
            varName_tmp = model_desc_sub{nVarMDLIdxs(jn), 1};
            
            indSpecVarIdx = find(ismember(trimSpec_in(:,1), varName_tmp)); % Has this specific variable name been individually specified?
            
            if isempty(indSpecVarIdx) %It has not been individually specified
                trimSpec_tmp = [trimSpec_tmp; {varName_tmp, nVarValue}]; % Add a row with the 'generic' value
            else
                trimSpec_tmp = [trimSpec_tmp; {varName_tmp, trimSpec_in{indSpecVarIdx,2}}]; % Add a row with the individually specified value
                
                rowsToBeDeleted = [rowsToBeDeleted; indSpecVarIdx]; % Make sure the individually specified row will be removed
            end
        end
                       
    end
    
    trimSpec_out(rowsToBeDeleted,:) = [];
    trimSpec_out = [trimSpec_out; trimSpec_tmp];

end
