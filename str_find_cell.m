function matching_cells = str_find_cell(cell_array, substring)
% Finds cells that contain the substring.
%
% Input:
%   - cell_array: Input cell array.
%   - substring: Substring to search for.
%
% Output:
%   - matching_cells: Logical array indicating cells containing the substring.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matching_cells = cellfun(@(x) ~isempty(strfind(x, substring)), ...
        cell_array, 'UniformOutput', true)';
end
