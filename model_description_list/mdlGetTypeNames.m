function type_names = mdlGetTypeNames(model_description,type)

switch type
    case 'state'
        if isempty(model_description.States.Continuous)
            cont = [];
        else
            cont = model_description.States.Continuous(:,1);
        end
        if isempty(model_description.States.Discrete)
            disc = [];
        else
            disc = model_description.States.Discrete(:,1);
        end
        type_names = [cont;disc];
    case 'input'
        type_names = model_description.MDL_Inputs;
    case 'output'
        type_names = model_description.MDL_Outputs;
    case 'deriv'
        type_names = [model_description.States.Continuous(:,1); ...
            model_description.States.Discrete(:,1)];
end

end

