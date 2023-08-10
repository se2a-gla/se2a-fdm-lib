function base_workspace_backup = baseWorkspaceBackup()
% baseWorkspaceBackup store current workspace in one variable
%   This function is useful if Simulink simulations are called from
%   functions using the sim function.
%   In case the Simulink model makes use of Simulink bus objects, the user
%   is forced to initialize the Simulink bus objects in the base workspace.
%   In this case, all other parameters of the Simulink model should be
%   located in the base workspace (not in the local function workspace).
%   Therefore, this function can be called at the beginning of the function
%   which calls the sim function. Then, all Simulink model parameters need
%   to be assigned to the base workspace using the assign function. Then,
%   the sim function should be called.
%   After extracting all Simulink model outputs, the base workspace can be
%   restored using the baseWorkspaceRestore function.
base_workspace_vars_list = evalin('base','whos');
base_workspace_backup = {};
for i = 1:length(base_workspace_vars_list)
    var_name = base_workspace_vars_list(i).name;
    base_workspace_backup{1,i} = var_name;
    base_workspace_backup{2,i} = evalin('base',var_name);
end
end