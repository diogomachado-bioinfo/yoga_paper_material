function hashes = words_to_hashinrecs(words, varargin)
% Words to Hash-inrec Conversion
%
% This function transforms words into hash-inrec values.
%
% Parameters:
%   words: Input words to be converted
%   varargin: Additional parameter-value pairs
%     - 'AMPLITUDE': Amplitude factor (default: 100000000000)
%     - 'MOD': Modulo factor (default: 1000000)
%
% Output:
%   hashes: Hash-inrec values obtained from the input words
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Handle optional input parameters
amplitude_factor = find_in_varargin(varargin, 'AMPLITUDE', 100000000000);
modulo_factor = find_in_varargin(varargin, 'MOD', 1000000);

% Check if the input words are in cell format
if iscell(words)
    % Convert cell array of words to a column vector
    words_vector = ensure_column_vector(words);
    
    % Apply string_to_recursive_indexing to each word, then compute hash-inrec values
    hashes = fix(mod(amplitude_factor * cellfun(@(x) string_to_recursive_indexing(x), words_vector, 'UniformOutput', true), modulo_factor));
else
    % Apply string_to_recursive_indexing to the input word, then compute hash-inrec value
    hashes = fix(mod(amplitude_factor * string_to_recursive_indexing(words), modulo_factor));
end

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

function recognition_results = string_to_recursive_indexing(str)
% String to Recursive Indexing Conversion
%
% This function converts a string into a DNA sequence, which is used for
% aligning sequences of characters.
%
% Parameters:
%   str: Input string to be converted
%
% Output:
%   recognition_results: Recognition results obtained using recursive indexing

% Get the number of rows in the input string matrix
num_rows = length(str(:, 1));

% Initialize an array to store recognition results
recognition_results = zeros(num_rows, 1);

% Process each row of the input string matrix in parallel
parfor i = 1:num_rows
    % Convert characters to binary, then to a matrix of 2 columns
    cols = vector_to_matrix(matrix_to_vector(de2bi(double(str(i, :)))), 2);
    
    % Convert binary digits to a DNA sequence, apply recursive indexing, and store the result
    recognition_results(i) = sequence_recursive_indexing(digits_to_dna((cols(:, 2) * 2 + cols(:, 1))'));
    
    % Display progress for every 10000 iterations
    if mod(i, 10000) == 0
        disp(i);
    end
end

end

function vector_result = matrix_to_vector(matrix_input)
% Matrix to Vector Conversion
%
% This function converts a matrix to a vector by concatenating its rows.
%
% Parameters:
%   matrix_input: Input matrix to be converted to a vector
%
% Output:
%   vector_result: Vector representation of the input matrix

% Transpose the input matrix to access its rows
matrix_transposed = matrix_input';

% Get the size of the transposed matrix
matrix_size = size(matrix_transposed);

% Initialize a vector to store the result
if ~iscell(matrix_input)
    % For numeric matrices, initialize with zeros
    result_vector = zeros(1, prod(matrix_size));
    
    % Assign matrix values to the vector, considering only non-zero elements
    result_vector(1:prod(matrix_size)) = matrix_transposed(matrix_transposed | ~matrix_transposed)';
else
    % For cell arrays, initialize with ones
    matrix_mask = ones(matrix_size(1), matrix_size(2));
    
    % Assign matrix values to the vector, considering only non-empty cells
    result_vector = matrix_transposed(matrix_mask == 1)';
end

% Return the resulting vector
vector_result = result_vector;

end

function matrix_result = vector_to_matrix(vector, width)
% Vector to Matrix Conversion
%
% This function creates a matrix with values from the input vector arranged by rows
% with the specified width.
%
% Parameters:
%   vector: Input vector to be converted to a matrix
%   width: Number of columns in the output matrix
%
% Output:
%   matrix_result: Matrix representation of the input vector

% Determine the number of elements in the input vector
num_elements = length(vector);

% Calculate the number of rows needed based on the specified width
num_rows = ceil(num_elements / width);

% Create an index matrix for accessing elements in the matrix
index_matrix = find(ones(width, num_rows));

% Initialize the output matrix
if iscell(vector)
    matrix = cell(width, num_rows);
else
    matrix = zeros(width, num_rows);
end

% Assign vector values to the matrix using the index matrix
matrix(index_matrix(1:num_elements)) = vector;

% Transpose the matrix to have the correct orientation
matrix_result = matrix';

end

function dna_sequence = digits_to_dna(digit_sequence)
% Digits to DNA Sequence Conversion
%
% This function converts numeric digits (0-3) to the corresponding DNA sequence.
%
% Parameters:
%   digit_sequence: Numeric digit sequence to be converted
%
% Output:
%   dna_sequence: DNA sequence represented by the input numeric digits

% Replace numeric digits with corresponding nucleotides
digit_sequence(digit_sequence == 0) = 'A';
digit_sequence(digit_sequence == 1) = 'C';
digit_sequence(digit_sequence == 2) = 'G';
digit_sequence(digit_sequence == 3) = 'T';

% Convert numeric digit sequence to character DNA sequence
dna_sequence = char(digit_sequence);

end

function recognition_result = sequence_recursive_indexing(seqs)
% Sequence Recursive Indexing Pattern Recognition
%
% This function performs pattern recognition on a DNA sequence using a recursive
% indexing approach.
%
% Parameters:
%   seqs: DNA sequences to be processed
%
% Output:
%   recognition_result: Recognition result obtained by the recursive indexing pattern recognition

% Convert DNA sequences to numeric digits
digit_sequences = dna_sequence_to_digits(seqs) + 1;

% Normalize digit values to the range [0, 1]
base_var = digit_sequences / 4;

% Apply recursive indexing pattern recognition
recognition_result = recursive_indexing_pattern_recognition(base_var, -1);

end

function digit_sequence = dna_sequence_to_digits(dna_sequence)
% DNA Sequence to Digits Conversion
%
% This function converts a DNA sequence to corresponding numeric digits (0-3).
%
% Parameters:
%   dna_sequence: Input DNA sequence to be converted
%
% Output:
%   digit_sequence: Numeric digit sequence representing the input DNA sequence

% Get the dimensions of the DNA sequence matrix
num_columns = length(dna_sequence(1, :));
num_rows = length(dna_sequence(:, 1));

% Identify positions of each nucleotide in the DNA sequence
cs = (upper(dna_sequence) == 'C');
gs = (upper(dna_sequence) == 'G');
ts = (upper(dna_sequence) == 'T');

% Note: Adenine (A) is not explicitly considered here, as it is represented by
% the absence of C, G, and T. The resulting digit will be 0 for adenine.

% Convert DNA sequence to numeric digits (0-3)
digit_sequence = zeros(num_rows, num_columns);
digit_sequence = digit_sequence + cs + 2 * gs + 3 * ts;

end

function recognition_result = recursive_indexing_pattern_recognition(weight_matrix, exponent)
% Recursive Indexing Pattern Recognition
% 
% This function performs pattern recognition using a recursive indexing approach
% based on the given weight matrix and exponent.
%
% Parameters:
%   weight_matrix: Weight matrix for pattern recognition
%   exponent: Exponent value for the recursive indexing function
%
% Reference:
% Souza, J. A.; "RECONHECIMENTO DE PADRÕES USANDO INDEXAÇÃO RECURSIVA";
% doctoral thesis; UFSC. pp 47.
%
% Example usage:
%   result = recursive_indexing_pattern_recognition(weight_matrix, 0.5);
%
% Output:
%   recognition_result: The final recognition result obtained by the algorithm.

% Initialize variables
num_columns = length(weight_matrix(1, :));
num_rows = length(weight_matrix(:, 1));
x = 5000000 * ones(num_rows, 1);
x_prime = 1;

% Recursive indexing process
for i = 1:num_columns
    x = sqrt((weight_matrix(:, num_columns - i + 1).^exponent) .* (0.0001 * tanh(x)));
end

% Final recognition result (averaging tanh values)
recognition_result = (tanh(x) + tanh(x_prime)) / 2;

end
