function PUBMED_STRUCT = parse_pubmed_table(pubmed_table, varargin)
% parse_pubmed_table - Parse PUBMED tables for text mining
%
%   PUBMED_STRUCT = parse_pubmed_table(pubmed_table, varargin) processes PUBMED tables for text mining purposes.
%
%   INPUT:
%       pubmed_table: Cell array of table filenames or preloaded tables.
%       varargin: Optional parameter-value pairs for customization.
%           - 'chcut': Character cut-off for text processing (default: 300).
%           - 'Tabin': Preloaded tables to use instead of loading from files (default: []).
%
%   OUTPUT:
%       PUBMED_STRUCT: Struct containing processed information and results.
%
%   EXAMPLE:
%       myPUBMED_STRUCT = parse_pubmed_table({'electrophoresis.tsv', 'mitochondrialgenome.tsv'});
%
%   NOTES:
%       - The function loads PUBMED tables, processes the title and abstract,
%         computes term frequencies, and generates various data structures for analysis.
%       - To run this function, it is necessary to download a required file
%         manually at: <https://bit.ly/4bJRazk>. The function expects this
%         file to be available at current folder or MATLAB path.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = varargin;
chrcut = find_in_varargin(V,'chcut',300); % Cut text according to character limit
tabin = find_in_varargin(V,'Tabin',[]);

process_num = 0;
process_num_max = 15;

% Load data set table
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Load data set table\n', process_num, process_num_max);
if isempty(tabin)
    all_text = [];
    for ii=1:length(pubmed_table)
        current_table = load_table_pubmed(pubmed_table{ii});
        all_text = [all_text; current_table];
    end
else
    all_text = tabin;
end

% Convert to uppercase and concatenate title and abstract
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Convert to uppercase and concatenate title and abstract\n', process_num, process_num_max);
processed_text = cellfun(@(X,Y) upper([X ' ' Y]), all_text(:,2), all_text(:,3), 'Un', 0);

% Remove null entries and filter by character length
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Remove null entries and filter by character length\n', process_num, process_num_max);
is_empty = (cellfun(@(x) isempty(x), all_text(:,3)));
text_lengths = (cellfun(@length, processed_text, 'Un',1));
processed_text = upper(processed_text(text_lengths > chrcut & ~is_empty));             
processed_table = all_text(text_lengths > chrcut & ~is_empty,:);
title_list = upper(processed_table(:,2));

% Keep unique title_list and corresponding text
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Keep unique title_list and corresponding text\n', process_num, process_num_max);
[~, unique_title_list] = unique(title_list);
processed_text = processed_text(unique_title_list);
title_list = title_list(unique_title_list);
processed_table = processed_table(unique_title_list,:);

% Process words from the concatenated text
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Process words\n', process_num, process_num_max);
ensured_words = ensure_column_vector(phrase_to_words(processed_text));
ensured_words = trim_all(only_letters_in_cell(ensured_words));
non_empty_cells = cellfun(@(x) ~isempty(x), ensured_words);
ensured_words = ensured_words(non_empty_cells);

% Obtain unique words and hash values
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Obtain unique words and hash values\n', process_num, process_num_max);
[unique_words, ~, unique_word_indices] = unique(ensured_words); 
mod_value = 1000000;
hash_in_texts = words_to_hashinrecs(unique_words,'MOD',mod_value);
hash_in_texts_all = hash_in_texts(unique_word_indices);
unique_words = ensure_column_vector(unique_words);
hash_in_texts = words_to_hashinrecs(unique_words,'MOD',mod_value);

% Identify unique hash values and initialize necessary data structures
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Identify unique hash values\n', process_num, process_num_max);
unique_hashes = unique(hash_in_texts);
num_unique_words = length(unique_hashes);
num_texts = length(processed_text);

hash_indices = cell(num_unique_words,1);
text_hashes = cell(num_texts,1);
word_list = cell(num_unique_words,1);
index_table = zeros(mod_value,1);
index_table(unique_hashes) = 1:num_unique_words;

% Selects when there are ambiguous hashes (highest occurrence)
process_num = process_num + 1;
reverse_str = '';
fprintf('\t[%02i/%02i] Selects when there are ambiguous hashes: ', process_num, process_num_max);
for i=1:num_unique_words
    occurrences = count_cell_occurrences(ensured_words((hash_in_texts_all==unique_hashes(i))));
    name_of_hash = occurrences.xcell{1};
    word_list{i} = name_of_hash;
    
    msg = sprintf('%i/%i', i, num_unique_words);
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n');

% Frequency of inrecs (word_freq) and words used for each hash (word_hash)
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Get inrecs frequency\n', process_num, process_num_max);
ctir = count_occurrences(hash_in_texts_all);
cttt = sum(ctir(:,2));
word_freq = ones(mod_value,1)/cttt;
word_freq(ctir(:,1)) = ctir(:,2)/cttt;
word_hash = cell(1000000,1);
word_hash(unique_hashes) = word_list;

% Normalize inrecs frequency
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Normalize inrecs frequency\n', process_num, process_num_max);
load RSMFRQ.mat % data source: https://www.kaggle.com/datasets/rtatman/english-word-frequency
normalized_freq = word_freq./RSMFRQ.FRQ;

% Generate indices of relationships between terms and documents
process_num = process_num + 1;
reverse_str = '';
fprintf('\t[%02i/%02i] Generate indices of relationships between terms and documents: ', process_num, process_num_max);
for i=1:num_texts
    words_in_text = trim_all(only_letters_in_cell(phrase_to_words(processed_text(i))));
    words_in_text = words_in_text(cellfun(@(x) ~isempty(x),words_in_text));
    hash_in_text = words_to_hashinrecs(words_in_text)';
    hash_indices(index_table(hash_in_text)) = cellfun(@(x) [x i], hash_indices(index_table(hash_in_text)), 'Un', 0);
    text_hashes{i} = hash_in_text;

    msg = sprintf('%i/%i', i, num_texts);
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n');

% Import base SWeeP
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Import base SWeeP\n', process_num, process_num_max);
load sweep-default-projection-matrix-1369.mat;
base_sweep = R1369;
num_sweep_rows = length(base_sweep(1,:));

% Sweep texts: DNA -> 3 frames AA
process_num = process_num + 1;
reverse_str = '';
fprintf('\t[%02i/%02i] Sweep texts: ', process_num, process_num_max);
sweep_result = zeros(num_texts, num_sweep_rows);
for i=1:1000:num_texts
    sweep_result(i:min(i+1000, num_texts), :) = text_to_w3f(...
        processed_text(i:min(i+1000, num_texts)), base_sweep);

    msg = sprintf('%i/%i', i, 1000);
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
end
msg = sprintf('%i/%i', 1000, 1000);
fprintf([reverse_str, msg]);
fprintf('\n');

% Process publication years
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Process publication years\n', process_num, process_num_max);
try
    publication_years = cellfun(@(x) num2str(char(x(1:4))), processed_table(:,5), 'Un', 0);
    publication_years = cellfun(@(x) str2num(x), publication_years, 'Un', 1);
catch
    publication_years = cell2mat(processed_table(:,5));
end

% Assign output struct
process_num = process_num + 1;
fprintf('\t[%02i/%02i] Assign output struct\n', process_num, process_num_max);
PUBMED_STRUCT.processed_table = processed_table;
PUBMED_STRUCT.processed_text = processed_text;
PUBMED_STRUCT.title_list = title_list;
PUBMED_STRUCT.text_lengths = text_lengths;
PUBMED_STRUCT.word_list = word_list;
PUBMED_STRUCT.index_table = index_table;
PUBMED_STRUCT.unique_hashes = unique_hashes;
PUBMED_STRUCT.word_hash = word_hash;
PUBMED_STRUCT.word_freq = word_freq;
PUBMED_STRUCT.hash_indices = hash_indices;
PUBMED_STRUCT.text_hashes = text_hashes;
PUBMED_STRUCT.num_texts = num_texts;
PUBMED_STRUCT.num_unique_words = num_unique_words;
PUBMED_STRUCT.publication_years = publication_years;
PUBMED_STRUCT.normalized_freq = normalized_freq;
PUBMED_STRUCT.base_sweep = base_sweep;
PUBMED_STRUCT.num_sweep_rows = num_sweep_rows;
PUBMED_STRUCT.sweep_result = sweep_result;
PUBMED_STRUCT.RSMFRQ = RSMFRQ;
PUBMED_STRUCT.all_text = all_text;
end

function modified_str = only_letters_in_cell(xstr, varargin)
% Only Letters in Cell
%
% This function replaces non-letter characters in xstr (string or cell array of strings)
% with spaces. If xstr is a cell array, the operation is applied to all elements.
%
% Parameters:
%   xstr: Input string or cell array of strings to be processed
%   varargin: Optional parameter indicating whether to include numbers (1 for true)
%
% Output:
%   modified_str: Processed string or cell array of strings with only letters (and optionally numbers)

% Initialize optional parameter
include_numbers = [];
if ~isempty(varargin)
    include_numbers = varargin{1};
end

% Check if xstr is a cell array
if iscell(xstr)
    % Apply only_letters to each element in the cell array
    modified_str = trim_all(cellfun(@(x) only_letters(x, include_numbers), xstr, 'UniformOutput', false));
else
    % Apply only_letters to the input string
    modified_str = only_letters(xstr, include_numbers);
end

end

function result_struct = count_cell_occurrences(input_cell_array)
% COUNT_CELL_OCCURRENCES - Count occurrences of unique elements in a cell array.
%   result_struct = COUNT_CELL_OCCURRENCES(input_cell_array)
%
%   INPUT:
%       - input_cell_array: Cell array of elements to count occurrences.
%
%   OUTPUT:
%       - result_struct: Structure containing occurrence information.

% Transform input cell array using ensure_column_vector function
input_cell_array = ensure_column_vector(input_cell_array);

% Find unique elements in the input cell array
u = unique(input_cell_array);
us = u;

% Get the number of unique elements
n = length(u);

% Initialize a matrix 'w' with the first column as indices and the second column with zeros
w = [(1:n)' zeros(n,1)];

% Loop through each unique element and count occurrences in the target cells
for i=1:n
    w(i,2) = sum(count_element_occurrences_in_target_cells(u(i),input_cell_array));
    us{i} = ['Option ' num2str(w(i,1)) ': ' u{i} ' - ' num2str(w(i,2)) ' times'];
end

% Sort the indices based on the second column of matrix 'w' in descending order
ids = sort_indices(w(:,2),'descend');

% Create a result struct with sorted information
result_struct.ocorr = w(ids,:);
result_struct.xcell = u(ids);
result_struct.strings = us(ids);

end

function occurrence_count = count_element_occurrences_in_target_cells(cell_array1, cell_array2)
% COUNT_OCCURRENCES_IN_CELLS - Count occurrences of elements from one cell array in another.
%   occurrence_count = COUNT_OCCURRENCES_IN_CELLS(cell_array1, cell_array2)
%
%   INPUTS:
%       - cell_array1: Cell array of elements to count occurrences.
%       - cell_array2: Cell array in which occurrences are counted.
%
%   OUTPUT:
%       - occurrence_count: Vector with the count of occurrences for each element in cell_array1.

% Use cellfun to count occurrences of each element in cell_array1 within cell_array2
occurrence_count = cellfun(@(x) sum(ismember(cell_array2, x)), cell_array1);
end

function W_matrix = text_to_w3f(text_cells, R)
% text_to_w3f - Converts text cells (TXT) into W matrix with Sweep projection in base R
% Originated from the conversion of 3 frames of DNA to amino acid sequences
%
% Parameters:
% - text_cells: Cell array of texts (TXT)
% - R: Projection base
%
% Output:
% - W_matrix: Resulting matrix

% Convert text cells to a structure with amino acid sequences
sequence_structure = cell2struct([cell(1, length(text_cells)); ...
    ensure_column_vector(cellfun(@(x) convert_dna_to_amino_acid_3frames(dnabits(x)), ...
    text_cells, 'UniformOutput', 0))'], {'Header', 'Sequence'});

% Generate W matrix using the Sweep projection in base R
W_matrix = generate_W160k_vectors(sequence_structure) * R;
end

function amino_acid_sequence = convert_dna_to_amino_acid_3frames(input_dna_sequence)
% convert_dna_to_amino_acid_3frames - Transforms DNA sequence into amino
% acid sequence in 3 frames separated by '****'
%
% Parameters:
% - input_dna_sequence: Input DNA sequence
%
% Output:
% - amino_acid_sequence: Resulting amino acid sequence

% Convert DNA sequence to amino acid sequence in 3 frames, separated by '****'
amino_acid_sequence = [nt2aa(input_dna_sequence(1:end), 'ACGTOnly', false) '****' ...
                      nt2aa(input_dna_sequence(2:end), 'ACGTOnly', false) '****' ...
                      nt2aa(input_dna_sequence(3:end), 'ACGTOnly', false)];

end

function W160k_vectors = generate_W160k_vectors(xfas, varargin)
% generate_W160k_vectors - Transform sequence structures into a matrix of 160k vectors
% Optionally, the matrix can be reordered according to the specified order (iord)
%
% Parameters:
% - xfas: Sequence structures
% - varargin: Optional parameters (e.g., 'Withpos', 'toSum', 'Ordem')
%
% Output:
% - W160k_vectors: Resulting matrix with 160k vectors

% Parse optional parameters
V = varargin;
withPos = find_in_varargin(V, 'Withpos', 0);
toSum = find_in_varargin(V, 'toSum', 0);
iord = find_in_varargin(V, 'Ordem', 1:160000);

% Generate prime numbers for indexing
prime_indices = generate_primes_up_to_n(1);
max_sequence_length = max(calculate_sequence_lengths(xfas));

if withPos
    prime_indices = generate_primes_up_to_n(max_sequence_length);
end

% Extract sequence structures
sequence_cell_array = struct2cell(xfas);
sequences = sequence_cell_array(2, :);

% Convert sequences to matrices and then to inline vectors
matrices = cellfun(@(x) matrix_to_inline(aa_sequence_to_matrix(...
    x, prime_indices, withPos, toSum)), sequences, 'UniformOutput', false);

% Concatenate matrices into a single matrix
W160k_vectors = cell2mat(matrices');
W160k_vectors = W160k_vectors(:, iord);

end

function [matrix_result, coordinates] = aa_sequence_to_matrix(x_sequence, Ps, with_pos, to_sum)
% x_sequence - amino acid sequence
% Ps - prime numbers for indexing (optional)
% with_pos - flag to include positional information
% to_sum - flag to sum occurrences instead of setting to 1

% Default sampling length (e.g., 1 for monopeptide, 2 for dipeptide, ...)
l = 2;

% Calculate the sampling space size
s = 20^l;

% Convert DNA sequence to a matrix
L5 = dna_sequence_to_matrix(x_sequence, l*2+1);

% Extract coordinates from amino acid numeric representations
coordinates = [amino_acid_to_numeric(L5(:,1:l)) amino_acid_to_numeric(L5(:,l+2:l*2+1))];

% Remove invalid coordinates
invalid_coordinates = ~prod(1 - double(coordinates <= 0), 2);
coordinates(invalid_coordinates, :) = [];

% Convert coordinates to indices
indices = coordinates_to_indices(coordinates, s);
indices(indices > (s^2)) = [];

% Initialize the matrix
matrix_result = zeros(s, s, 'double');

if with_pos
    % Include positional information using prime numbers
    if nargin < 2 || isempty(Ps) || length(Ps) < length(indices)
        Ps = generate_primes_up_to_n(length(indices));
    end
    
    % Unify identical indices using prime numbers and calculate values
    [values, indices] = unifysameinds(indices, Ps(1:length(indices)), @(x) prod(x.^(1/length(x))));
    matrix_result(indices) = values;
else
    if to_sum
        % Sum occurrences for identical indices
        [values, indices] = unifysameinds(indices, ones(1, length(indices)), @(x) sum(x));
        matrix_result(indices) = values;  % You can also use log(vals+1) here
    else
        % Set occurrences to 1 without summing
        matrix_result(indices) = 1;
    end
end

end

function prime_numbers = generate_primes_up_to_n(n)
% Returns the first n prime numbers (2:n)

% Initialize variables
m = n;
n_max = 100;
primes_found = find(isprime(1:n_max));

% Find primes until the desired count is reached
while length(primes_found) < m
    n_max = n_max * 2;
    primes_found = find(isprime(1:n_max));
end    

% Return the result
prime_numbers = primes_found(1:m);

end

function sequence_lengths = calculate_sequence_lengths(structure_array)
% Calculates the lengths of sequences in an array of structures.
% Each structure in the array should have a 'Sequence' field.

% Get the number of structures in the array
num_structures = length(structure_array);

% Initialize a vector to store the lengths of sequences
sequence_lengths = zeros(num_structures, 1);

% Loop through each structure and calculate the length of its sequence
for i = 1:num_structures
    sequence_lengths(i) = length(structure_array(i).Sequence);
end

end

function pattern_matrix = dna_sequence_to_matrix(dna_segment, pattern_length)
% Creates a matrix 'pattern_matrix' with patterns (examples) per row for a
% given segment of DNA

% Calculate the length of the DNA segment
n = length(dna_segment);

% Generate column indices for the resulting matrix
column_indices = uint32(1:(n-pattern_length+1))';

% Create indices matrix for extracting patterns
indices = repmat(column_indices, 1, pattern_length) + ...
    repmat(uint32(1:pattern_length), (n-pattern_length+1), 1) - 1;

% Extract patterns from the DNA segment using indices
pattern_matrix = dna_segment(indices);

% Transpose the matrix if the pattern length is 1
if pattern_length == 1
    pattern_matrix = pattern_matrix';
end

end

function numeric_representation = amino_acid_to_numeric(input_sequence)
% Convert amino acid sequence 'input_sequence' to numeric representation.
% Each amino acid is assigned a unique numeric value based on its order.

% Get the size of the input sequence
[num_sequences, sequence_length] = size(input_sequence);

% Convert amino acids to integers and subtract 1
values = double(aa2int(upper(input_sequence))) - 1;

% Create a matrix of powers for subsequent calculations
powers = repmat(double(0:(sequence_length - 1)), num_sequences, 1);

% Calculate the numeric representation using base 20 conversion
numeric_representation = sum((repmat(20, num_sequences, sequence_length).^powers) .* values, 2) + 1;

% Identify invalid amino acids and set corresponding numeric values to -1
invalid_amino_acids = ~prod(1 - double(values < 0 | values > 19), 2);
numeric_representation(invalid_amino_acids) = -1;
    
end

function matrix_indices = coordinates_to_indices(ij_coordinates, total_columns)
% Generate matrix indices for a scan pattern given pairs of (i, j) coordinates and the total number of columns 'total_columns'.
% Each pair in 'ij_coordinates' corresponds to (row i, column j), and the resulting indices are arranged with 'total_columns' elements per column.

matrix_indices = (ij_coordinates(:, 2) - 1) * total_columns + ij_coordinates(:, 1);
    
end

function output_vector = matrix_to_inline(input_matrix)
% Linearize a matrix 'input_matrix' by columns.
% The resulting vector 'output_vector' contains elements from 'input_matrix' arranged by columns.

output_vector = input_matrix(input_matrix == input_matrix)';

end
