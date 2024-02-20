function sorted_indices = sort_indices(vect, varargin)
% SORT_INDICES - Return the indices that would sort a vector in ascending
% or descending order.
%   sorted_indices = SORT_INDICES(vect, varargin)
%
%   INPUTS:
%       - vect: Input vector to be sorted.
%       - varargin: Optional argument specifying sorting mode ('ascend' or 'descend').
%
%   OUTPUTS:
%       - sorted_indices: Indices that would sort the input vector.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize sorting mode if provided
if ~isempty(varargin)
    sorting_mode = varargin{1};
end

% Sort the vector and get the indices
if isempty(varargin)
    [~, sorted_indices] = sort(vect);
else
    [~, sorted_indices] = sort(vect, sorting_mode);
end

end
