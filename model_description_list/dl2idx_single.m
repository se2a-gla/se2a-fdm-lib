function [varIdxs] = dl2idx_single(list, varNames, col)
     
%Function to convert a single variable name from the model description list into its respective index/indices. 
% NB: variable name can be non-unique - all matching indices will be returned in varIdxs.
%
% Inputs: -list:     Lowest substructure of model description list struct, e.g modelDescription.AC_Outputs
%         -varNames: Matrix containing strings of names corresponding to the chosen column in the description list
%         -col:      Desired column index in description list.

[nbRows, nbCols] = size(varNames);

if iscell(varNames) && (nbRows > 1 || nbCols > 1)
    error('Too many variable names entered!')
end

if iscell(varNames)
    nChar = length(varNames{1});
else
    nChar = length(varNames);
end

varIdxs = find(strncmp(varNames, list(:,col), nChar));
end