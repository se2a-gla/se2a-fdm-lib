function trim_point = tpSetMdlNames(trim_point,type,sub_names,varargin)

switch type
    case 'state'
        if isempty(varargin)
            num_cont = length(trim_point.mdl_names.States.Continuous);
        else
            num_cont = varargin{1};
        end
        trim_point.mdl_names.States.Continuous = sub_names(1:num_cont);
        trim_point.mdl_names.States.Discrete = sub_names(num_cont+1:end);
    case 'input'
        trim_point.mdl_names.MDL_Inputs = sub_names;
    case 'output'
        trim_point.mdl_names.MDL_Outputs = sub_names;
    case 'deriv'
        trim_point.mdl_names.Deriv = sub_names;
end

end

