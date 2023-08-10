function StateMDL_final = mdlConstructState(statePathList, StateMDL_red, numStates, varargin)

    if isempty(varargin)
        idx_begin = 1;
    else
        idx_begin = varargin{1};
    end
    
    StateMDL_final = cell(numStates, max(1,size(StateMDL_red,2))); % Preallocate final MDL
    stateCount = 0;
    unmatchedStates = [];

    for in = 1:size(StateMDL_red, 1)

        idx_tmp = find(strcmp(StateMDL_red{in,5}, statePathList));
        stateCount = stateCount + size(idx_tmp,1);

        if isempty(idx_tmp) % Do none of the states match?
%             warning('Expected states not found in model');
%             unmatchedStates = [unmatchedStates; StateMDL_red{in,1}];
            unmatchedStates = [unmatchedStates; StateMDL_red(in,1)];
            continue;
        else
            col_tmp = cell(size(idx_tmp,1),4);
            for kn = 1:4 % Cycle through columns 1-4
                if (size(idx_tmp,1)>1) && (size(StateMDL_red{in,kn},1)==1) %Is only 1 entry defined for multiple states?
                    if (kn == 1)      % If in columns 1 or 3, number the copies.
                        col_tmp(:,kn) = strcat(StateMDL_red(in,kn), '__', num2str([1:size(idx_tmp,1)]','%-u'));
                    elseif (kn == 3)
                        col_tmp(:,kn) = strcat(StateMDL_red(in,kn), '_', num2str([1:size(idx_tmp,1)]','%-u')); %Only 1 underscore in col. 3 for better formatting
                    else              % Otherwise, just copy them.
                        col_tmp(:,kn) = StateMDL_red(in,kn);
                    end

                elseif size(idx_tmp,1)==size(StateMDL_red{in,kn},1) % Do the numbers of specified states match existing states?
                    if ~iscell(StateMDL_red{in,kn}) %Is there only one entry? (i.e. not a cell array)
                        col_tmp(:,kn) = StateMDL_red(in,kn);
                    else
                        col_tmp(:,kn) = StateMDL_red{in,kn};
                    end

                elseif size(StateMDL_red{in,kn},1) < 1 %Is this entry empty?
                    % Do nothing.

                else % In this case, more than one state has been specified, 
                     % but it does not match the number of states with the same path reported by the model.
                    error('Mismatch in specified and reported number of states');
                end
            end
        end

        for jn  = 1:size(idx_tmp,1) % Assign values to final MDL
            StateMDL_final(idx_tmp(jn), 1:4) = col_tmp(jn,1:4);
            StateMDL_final{idx_tmp(jn), 5} = StateMDL_red{in,5};
        end

    end
    
    if ~isempty(unmatchedStates) % Were any states defined in the MDL not found in the model?
        warning('The following defined states were not found in the model:');
        for jn=1:length(unmatchedStates)
            disp(unmatchedStates{jn});
        end
    end
    
    if stateCount < numStates % Did all states loaded above cover all states declared by the model?
        
        idx_missingStates = find(cellfun(@isempty, StateMDL_final(:,1)));
        num_missingStates = idx_missingStates + idx_begin - 1;
        
        generic_stateNames = generate_numbered_str_array('x__',num_missingStates);

        warning('The following states are not accounted for in Model Description List and have been assigned generic names:');
        disp('States paths:');
        disp(statePathList(idx_missingStates));
        disp('Generic names:');
        disp(generic_stateNames);
                
        for in = 1:length(idx_missingStates)
            StateMDL_final{idx_missingStates(in),1} = generic_stateNames{in,:};
            StateMDL_final(idx_missingStates(in),5) = statePathList(idx_missingStates(in),1);
        end
    end
end