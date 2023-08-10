function [mdl_nbd] = name_MultipleAileronPairs(n_ailerons, mdl_row)
% Name multiple pairs of ailerons, from outermost to innermost. 
% Inputs:
% - n_ailerons: Number of aileron pairs
% - mdl_row: (Optional) Prototypical row of the model description list which 
%               is filled in and extended as necessary. 
%               If only the first column is provided e.g. {'da_sym'}, mdl_ndb only contains
%               the first column of the mdl section.
%               If nothing is provided, mdl_row simply returns the
%               suffixes for the first column.

if nargin<2
    mdl_row = {''};
end

n_col = size(mdl_row, 2);

switch n_ailerons
    case 1 % Only 1 pair, no naming
        mdl_nbd = mdl_row;
    case 2 % Inner, Outer
        mdl_nbd = repmat(mdl_row, n_ailerons, 1);
        suffix_set = {'__out', '^[o]', ', outer'
                   '__in', '^[i]', ', inner'; 
                   };

    case 3 % Inner, Middle, Outer
        mdl_nbd = repmat(mdl_row, n_ailerons, 1);
        suffix_set = {'__out', '^[o]', ', outer'; 
                   '__mid', '^[m]', ', middle'; 
                   '__in', '^[i]', ', inner'};

    otherwise % Numbered (from outermost to innermost)
        mdl_nbd = repmat(mdl_row, n_ailerons, 1);
        suffix_set = cell(n_ailerons, 3);%[];
        for ii = 1:n_ailerons
            suffix_set{ii,1} = ['__', num2str(ii)];
            suffix_set{ii,2} = ['^[', num2str(ii),']'];
            suffix_set{ii,3} = [', n.', num2str(ii)];
        end
end

for ii = 1:n_ailerons
    mdl_nbd{ii,1} = [mdl_row{1,1}, suffix_set{ii,1}];
    if n_col > 1
        mdl_nbd{ii,3} = [mdl_row{1,3}, suffix_set{ii,2}];
        mdl_nbd{ii,4} = [mdl_row{1,4}, suffix_set{ii,3}];
    end
end

end


