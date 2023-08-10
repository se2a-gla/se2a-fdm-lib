function mdl_names = tpGetMdlNames(trim_point,type)

switch type
    case 'state'
        mdl_names = [ trim_point.mdl_names.States.Continuous; ...
                    trim_point.mdl_names.States.Discrete ];
    case 'input'
        mdl_names = trim_point.mdl_names.MDL_Inputs;
    case 'output'
        mdl_names = trim_point.mdl_names.MDL_Outputs;
    case 'deriv'
        mdl_names = trim_point.mdl_names.Deriv;
end

end

