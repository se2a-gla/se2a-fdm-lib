function trim_point = tpRmVal(trim_point,type,val_names)
% tpRmVal maniuplate trim point: remove specified value(s)

if ~isa(val_names,'cell')
    val_names = {val_names};
end

sub_names = mdlGetTypeNames(trim_point.mdl_names,type);
sub_array = tpGetType(trim_point,type);

for i = 1:length(val_names)
    if strcmp(val_names{i}(end-1:end),'__')
        num_removed = 0;
        while true
            is_removed = false;
            for j = 1:length(sub_names)
                if contains(sub_names{j},val_names{i})
                    trim_point = tpRmVal(trim_point,type,sub_names(j));
                    sub_names = mdlGetTypeNames(trim_point.mdl_names,type);
                    is_removed = true;
                    num_removed = num_removed + 1;
                    break;
                end
            end
            if ~is_removed
                % disp(['Removed ',num2str(num_removed),' values for ',val_names{i},'.']);
                break;
            end
        end
    else
        idx = dl2idx(sub_names, val_names(i));
        if idx~=0
            sub_names(idx) = [];
            sub_array(idx) = [];
            num_cont = length(trim_point.mdl_names.States.Continuous);
            if idx <= num_cont
                num_cont = num_cont - 1;
            end
            trim_point = tpSetType(trim_point,type,sub_array);
            trim_point = tpSetMdlNames(trim_point,type,sub_names,num_cont);
        end
    end
end

end

