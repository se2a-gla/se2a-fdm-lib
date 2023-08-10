function sub_array = tpGetType(trim_point,type)

switch type
    case 'state'
        sub_array = trim_point.State;
    case 'input'
        sub_array = trim_point.Input;
    case 'output'
        sub_array = trim_point.Output;
    case 'deriv'
        sub_array = trim_point.Deriv;
end

end

