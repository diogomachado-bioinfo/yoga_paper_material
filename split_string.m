function cell_elements = split_string(s, varargin)
% Split String
%
% This function returns a cell array with elements of the input string separated
% by space or a specified character in varargin.
%
% Parameters:
%   s: Input string to be split
%   varargin: Optional parameter specifying the character for splitting
%
% Output:
%   cell_elements: Cell array with elements of the input string
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the default splitting character to space (32)
split_character = ' ';

% Check if varargin is not empty, if yes, use the specified character
if ~isempty(varargin)
    split_character = varargin{1};
end

% Add a space at the end to ensure the last element is included in the split
string_with_spaces = [eliminate_repetitions(s) split_character];

% Find indices where the string should be split
split_indices = [1 find(string_with_spaces == split_character | string_with_spaces == 10) + 1];

% Initialize cell array to store the split elements
cell_elements = cell(1, length(split_indices) - 1);

% Loop through the split indices and extract each element
for i = 1:length(split_indices) - 1
    cell_elements{i} = string_with_spaces(split_indices(i):split_indices(i + 1) - 2);
end

end
