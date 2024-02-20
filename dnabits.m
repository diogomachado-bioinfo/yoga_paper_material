function dna_sequence = dnabits(input_string)
% Convert a string 'input_string' to a DNA-like sequence.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert the input string to a matrix of bits
bits_matrix = byte_to_bits(input_string);

% Convert the bits matrix to a vector
bits_vector = matrix_to_vector(bits_matrix);

% Reshape the vector into a matrix with 2 columns
reshaped_matrix = vector_to_matrix(bits_vector, 2);

% Combine the two columns into a single numeric vector
numeric_sequence = reshaped_matrix(:, 2) * 2 + reshaped_matrix(:, 1);

% Convert the numeric vector to a DNA-like sequence
dna_sequence = numeric_to_dna(numeric_sequence');

end

function output_bits = byte_to_bits(input_bytes)
% Convert an array of bytes 'input_bytes' to a matrix of bits.
% Each byte is represented by a row in the output matrix 'output_bits',
% with the individual bits as columns.

% Convert the input bytes to double
input_bytes = double(input_bytes);

% Get the length of the input vector
num_bytes = length(input_bytes);

% Initialize the output matrix with zeros
output_bits = zeros(num_bytes, 8);

% Create a temporary variable to store the original byte values
temp_bytes = input_bytes;

% Loop through each bit position (from right to left)
for i = 1:8
    % Identify bits at the current position
    bits_at_position = ~(fix(temp_bytes / 2) == (temp_bytes / 2));

    % Update the output matrix with the identified bits for the current position
    output_bits(bits_at_position, i) = 1;

    % Update the temporary variable by shifting right
    temp_bytes = fix(temp_bytes / 2);
end

end

function output_vector = matrix_to_vector(input_matrix)
% Convert a matrix 'input_matrix' into a vector, arranging elements by rows.

% Transpose the input matrix to access elements by columns
transposed_matrix = input_matrix';

% Get the size of the transposed matrix
matrix_size = size(transposed_matrix);

% Check if the matrix is not a cell array
if ~iscell(input_matrix)
    % Create a vector of zeros with the size equal to the product of matrix dimensions
    vector_elements = zeros(1, prod(matrix_size));

    % Assign non-zero elements from the transposed matrix to the vector
    vector_elements(1:prod(matrix_size)) = transposed_matrix(transposed_matrix | ~transposed_matrix)';

else
    % For cell arrays, create a matrix of ones with the same size as the transposed matrix
    ones_matrix = ones(matrix_size(1), matrix_size(2));

    % Extract non-zero elements from the transposed matrix and arrange them in a vector
    vector_elements = transposed_matrix(ones_matrix == 1)';

end

% Assign the resulting vector to the output variable
output_vector = vector_elements;

end

function result_matrix = vector_to_matrix(vect, larg)
% Convert a vector 'vect' into a matrix with 'larg' columns.
% The resulting matrix 'result_matrix' has values from 'vect' arranged by row,
% with the last row padded with zeros if necessary.

% Calculate the number of elements in the vector
num_elements = length(vect);

% Calculate the number of rows needed in the matrix
num_rows = ceil(num_elements / larg);

% Create indices to fill the matrix by column
indices = find(ones(larg, num_rows));

% Initialize the matrix with zeros or cells based on the input type
if iscell(vect)
    result_matrix = cell(larg, num_rows);
else
    result_matrix = zeros(larg, num_rows);
end

% Assign values from the vector to the matrix using the calculated indices
result_matrix(indices(1:num_elements)) = vect;

% Transpose the matrix to have values arranged by row
result_matrix = result_matrix';

end

function dna_sequence = numeric_to_dna(numeric_sequence)
% Convert a numeric vector representing DNA bases to a DNA sequence string.
% The input vector 'numeric_sequence' represents DNA bases (0=A, 1=C, 2=G, 3=T).
% The output 'dna_sequence' is a string containing the corresponding DNA sequence.

% Replace numeric values with corresponding DNA bases
numeric_sequence(numeric_sequence == 0) = 'A';
numeric_sequence(numeric_sequence == 1) = 'C';
numeric_sequence(numeric_sequence == 2) = 'G';
numeric_sequence(numeric_sequence == 3) = 'T';

% Convert the numeric sequence to a character string representing the DNA sequence
dna_sequence = char(numeric_sequence);

end
