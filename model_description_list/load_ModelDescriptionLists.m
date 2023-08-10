function [model_description] = load_ModelDescriptionLists(modelName)
%   Load model description lists.
%   Prefixes:
%   - MDL: Model

% load_system('sim_flexible_unsteady_dev');
load_system(modelName);

x0_off_on	= get_param(modelName,'LoadInitialState');
u_off_on	= get_param(modelName,'LoadExternalInput');

set_param(modelName,'LoadInitialState','off');
set_param(modelName,'LoadExternalInput','off');

model_description = [];

model_description.States = load_StatesMDL(modelName);

model_description.MDL_Inputs = load_InputMDL(modelName);

model_description.MDL_Outputs = load_OutputMDL(modelName, model_description.States);

set_param(modelName,'LoadInitialState',x0_off_on);
set_param(modelName,'LoadExternalInput',u_off_on);

close_system(modelName,0);

end

