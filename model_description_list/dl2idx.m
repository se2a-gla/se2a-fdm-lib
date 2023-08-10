function [varIdxs] = dl2idx(list, varNames, col, varargin)

%Function to convert a cell matrix of names from the model description list into their respective indexes.
% NB: this function assumes all submitted names are unique in the model description list.
%
% Inputs: -list:     Lowest substructure of model description list struct, e.g model_description.MDL_Outputs
%         -varNames: Cell matrix containing strings of names corresponding to the chosen column in the description list
%         -col:      (Optional)Desired column index in description list. If
%                    unspecified, column 1 is set by default.


if nargin  < 3
    col = 1;
end
if isempty(varargin)
    is_warn = true;
else
    is_warn = varargin{1};
end

[nbRows, nbCols] = size(varNames);
varIdxs = zeros(nbRows, nbCols);
for in = 1:nbRows
    for jn = 1:nbCols
        varName_tmp = varNames(in,jn);
        if iscell(varName_tmp)
            varName_tmp = char(varName_tmp);
        end
        
        match_tmp = find(strcmp(varName_tmp, list(:,col)));
        if size(match_tmp,1) == 1
            varIdxs(in,jn) = match_tmp;
        elseif isempty(match_tmp) && is_warn
            warning(['Requested MDL term "', varName_tmp, '" was not found.'])
        elseif is_warn
            warning(['Requested MDL term "',  varName_tmp,'" is not unique!']);
        end
    end
end
        