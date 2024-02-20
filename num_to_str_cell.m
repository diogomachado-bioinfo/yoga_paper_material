function string_cell = num_to_str_cell(numeric_vector)
% Numeric to String Cell
%
% This function creates a cell array of strings corresponding to the numbers in the input vector.
%
% Parameters:
%   numeric_vector: Numeric vector to be converted
%
% Output:
%   string_cell: Cell array of strings representing the numbers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure the input vector is a column vector
column_vector = ensure_column_vector(numeric_vector);

% Convert the numeric vector to a cell array of strings
string_cell = matrix_to_cell_lines(num2str(column_vector));
end

function column_vector = ensure_column_vector(input_vector)
% Ensure Column Vector
%
% This function ensures that the input vector is returned as a column vector.
%
% Parameters:
%   input_vector: Input vector to be processed
%
% Output:
%   column_vector: Column vector representation of the input

% Get the size of the input vector
vector_size = size(input_vector);

% Initialize the output as the input vector
column_vector = input_vector;

% Check if the vector is a row vector, then transpose it to make it a column vector
if vector_size(2) > vector_size(1)
    column_vector = input_vector';
end

end

function cell_lines = matrix_to_cell_lines(mat)
% Matrix to Cell Lines
%
% This function generates a cell array with each row of the matrix as a separate cell.
%
% Parameters:
%   mat: Input matrix to be converted
%
% Output:
%   cell_lines: Cell array with each row of the matrix as a separate cell

% Create a cell array with each row of the matrix as a separate cell
cell_lines = mat2cell(mat, ones(1, size(mat, 1)), size(mat, 2));
end
