function mdl_names = mdlGetMdlNames(model_description)

if ~isempty(model_description.States.Continuous)
    mdl_names.States.Continuous = model_description.States.Continuous(:,1);
else
    mdl_names.States.Continuous = [];
end
if ~isempty(model_description.States.Discrete)
    mdl_names.States.Discrete = model_description.States.Discrete(:,1);
else
    mdl_names.States.Discrete = [];
end

if ~isempty(model_description.MDL_Inputs)
    mdl_names.MDL_Inputs = model_description.MDL_Inputs(:,1);
else
    mdl_names.MDL_Inputs = [];
end

if ~isempty(model_description.MDL_Outputs)
    mdl_names.MDL_Outputs = model_description.MDL_Outputs(:,1);
else
    mdl_names.MDL_outputs = [];
end

end

