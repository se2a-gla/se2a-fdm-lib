function trim_point = tpSetType(trim_point,type,sub_array)

switch type
    case 'state'
        trim_point.State = sub_array;
    case 'input'
        trim_point.Input = sub_array;
    case 'output'
        trim_point.Output = sub_array;
    case 'deriv'
        trim_point.Deriv = sub_array;
end

end

