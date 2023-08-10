function [list_out] = generate_numbered_str_array(str_in,idx_set)
%GENERATE_NUMBERED_STR_ARRAY Accepts a leading string and a set of integer
%digits, and appends the digits to the string without extra spaces. 
%Outputs the list the the form of a column cell array.
idx_set = idx_set(:);

idx_set_cell = mat2cell(idx_set,ones(size(idx_set)));

list_out = strcat(str_in, cellfun(@num2str,idx_set_cell,'UniformOutput', false));

end

