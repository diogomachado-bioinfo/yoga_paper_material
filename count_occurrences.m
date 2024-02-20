function [occurrence_matrix] = count_occurrences(vector)
% COUNT_OCCURRENCES - Count the number of occurrences for each unique element in a vector.
%   [occurrence_matrix, display_list] = COUNT_OCCURRENCES(vector)
%
%   INPUTS:
%       - vector: Input vector.
%
%   OUTPUTS:
%       - occurrence_matrix: Matrix containing unique elements and their corresponding occurrences.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Count occurrences using count_unique function
[a, b] = count_unique(vector);
occurrence_matrix = [a, b];

% Sort the occurrence matrix in descending order based on occurrences
occurrence_matrix = occurrence_matrix(sort_indices(occurrence_matrix(:, 2), ...
    'descend'), :);

end

function [uniques,numUnique] = count_unique(x,option)
%COUNT_UNIQUE  Determines unique values, and counts occurrences
%   [uniques,numUnique] = count_unique(x)
%
%   This function determines unique values of an array, and also counts the
%   number of instances of those values.
%
%   This uses the MATLAB builtin function accumarray, and is faster than
%   MATLAB's unique function for intermediate to large sizes of arrays for integer values.  
%   Unlike 'unique' it cannot be used to determine if rows are unique or 
%   operate on cell arrays.
%
%   If float values are passed, it uses MATLAB's logic builtin unique function to
%   determine unique values, and then to count instances.
%
%   Descriptions of Input Variables:
%   x:  Input vector or matrix, N-D.  Must be a type acceptable to
%       accumarray, numeric, logical, char, scalar, or cell array of
%       strings.
%   option: Acceptable values currently only 'float'.  If 'float' is
%           specified, the input x vector will be treated as containing
%           decimal values, regardless of whether it is a float array type.
%
%   Descriptions of Output Variables:
%   uniques:    sorted unique values
%   numUnique:  number of instances of each unique value
%
%   Example(s):
%   >> [uniques] = count_unique(largeArray);
%   >> [uniques,numUnique] = count_unique(largeArray);
%
%   See also: unique, accumarray

% Author: Anthony Kendall
% Contact: anthony.kendall@gmail.com
% Created: 2009-03-17

testFloat = false;
if nargin == 2 && strcmpi(option,'float')
    testFloat = true;
end

nOut = nargout;
if testFloat
    if nOut < 2
        [uniques] = float_cell_unique(x,nOut);
    else
        [uniques,numUnique] = float_cell_unique(x,nOut);
    end
else
    try %this will fail if the array is float or cell
        if nOut < 2
            [uniques] = int_log_unique(x,nOut);
        else
            [uniques,numUnique] = int_log_unique(x,nOut);
        end
    catch %default to standard approach
        if nOut < 2
            [uniques] = float_cell_unique(x,nOut);
        else
            [uniques,numUnique] = float_cell_unique(x,nOut);
        end
    end
end

end

function [unique_values, num_unique_occurrences] = int_log_unique(input_data, n_out)
% INT_LOG_UNIQUE - Extract unique values and their occurrences from an array of integers.
%   [unique_values, num_unique_occurrences] = INT_LOG_UNIQUE(input_data, n_out)
%
%   INPUTS:
%       - input_data: Input array of integers.
%       - n_out: Number of outputs, set to 2 to also get the count of occurrences.
%
%   OUTPUTS:
%       - unique_values: Array containing unique integer values from the input.
%       - num_unique_occurrences: (Optional) Number of occurrences for each unique integer value.

% First, determine the offset for negative values
min_value = min(input_data(:));

% Check if accumarray is appropriate for this function
max_index = max(input_data(:)) - min_value + 1;
if max_index / numel(input_data) > 1000
    error('Accumarray is inefficient for arrays when index values are significantly larger than the number of elements.');
end

% Now, offset to get the index
index = input_data(:) - min_value + 1;

% Count the occurrences of each index value
num_unique_occurrences = accumarray(index, 1);

% Get the values which occur more than once
unique_indices = (1:length(num_unique_occurrences))';
unique_values = unique_indices(num_unique_occurrences > 0) + min_value - 1;

if n_out == 2
    % Trim the num_unique_occurrences array
    num_unique_occurrences = num_unique_occurrences(num_unique_occurrences > 0);
end

end

function [unique_values, num_unique_occurrences] = float_cell_unique(input_data, n_out)
% FLOAT_CELL_UNIQUE - Extract unique values and, optionally, the count of duplicate occurrences.
%   [unique_values, num_unique_occurrences] = FLOAT_CELL_UNIQUE(input_data, n_out)
%
%   INPUTS:
%       - input_data: Input array, can be either numeric or a cell array.
%       - n_out: Number of outputs, set to 2 to also get the count of duplicate values.
%
%   OUTPUTS:
%       - unique_values: Array containing unique values from the input.
%       - num_unique_occurrences: (Optional) Number of occurrences of each unique value.

% Check if the input is a cell array or a numeric array
if ~iscell(input_data)
    % For a numeric array, first sort the input vector
    sorted_data = sort(input_data(:));
    num_of_elements = numel(sorted_data);

    % Check if the array type needs to be converted to double
    original_class = class(sorted_data);
    is_double = isa(original_class, 'double');

    % Convert to double if necessary
    if ~is_double
        sorted_data = double(sorted_data);
    end

    % Check for NaNs or Infs; sort returns these at the beginning or end
    if isnan(sorted_data(1)) || isinf(sorted_data(1)) || ...
            isnan(sorted_data(num_of_elements)) || ...
            isinf(sorted_data(num_of_elements))
        % Check if the array contains NaNs or Infs
        nan_indices = isnan(sorted_data);
        inf_indices = isinf(sorted_data);
        indices_to_remove = nan_indices | inf_indices;

        % Remove NaNs and Infs from the array
        sorted_data = sorted_data(~indices_to_remove);
    end

    % Determine break locations of unique values
    unique_locations = [true; diff(sorted_data) ~= 0];
else
    is_double = true; % Just to avoid conversion at the end

    % For a cell array, sort the rows
    sorted_data = sort(input_data(:));

    % Determine unique location values
    unique_locations = [true; ~strcmp(sorted_data(1:end-1), sorted_data(2:end)) ~= 0];
end

% Determine the unique values
unique_values = sorted_data(unique_locations);

% Convert back to the original data type if it was not double
if ~is_double
    sorted_data = feval(original_class, sorted_data);
end

% Count the number of duplicate values if requested
if n_out == 2
    num_unique_occurrences = diff([find(unique_locations); length(sorted_data) + 1]);
end

end
