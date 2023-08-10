function trim_point = tpAddVal(trim_point,type2,val_names,vals)
% tpAddVal maniuplate trim point: add value at the end

% type2 is the type or 'continuous state' or 'discrete state' in case of
% states.


if ~isa(val_names,'cell')
    val_names = {val_names};
end

if strcmp(type2(end-4:end),'state')
    type = 'state';
else
    type = type2;
end

sub_names = mdlGetTypeNames(trim_point.mdl_names,type);
sub_array = tpGetType(trim_point,type);

for i = 1:length(val_names)
    if strcmp(val_names{i}(end-1:end),'__')
        for j = 1:length(vals)
            trim_point = tpAddVal(trim_point,type,{[val_names{i},num2str(j)]},vals(j));
        end
    else
        is_free = true;
        for j = 1:length(sub_names)
            if strcmp(sub_names{j},val_names{i})
                warning(['Variable ',val_names{i},' already exists.']);
                is_free = false;
                break;
            end
        end
        if ~is_free
            continue;
        end
        num_cont = length(trim_point.mdl_names.States.Continuous);
        if strcmp(type2,'continuous state')
            sub_names = [sub_names(1:num_cont);val_names(i);sub_names(num_cont+1:end)];
            sub_array = [sub_array(1:num_cont);vals(i);sub_array(num_cont+1:end)];
            num_cont = num_cont + 1;
        else
            sub_names(end+1) = val_names(i);
            sub_array(end+1) = vals(i);
        end
        trim_point = tpSetType(trim_point,type,sub_array);
        trim_point = tpSetMdlNames(trim_point,type,sub_names,num_cont);
    end
end

end

