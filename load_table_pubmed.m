function pubmed_table = load_table_pubmed(file_path)
% Loads a table downloaded from Pubmed with four fields:
% Column 1: PMID
% Column 2: Title
% Column 3: Abstract
% Column 4: Authors
%
% Inputs:
%   - file_path: The path to the file to be loaded.
%
% Output:
%   - pubmed_table: A cell array with four columns for each record.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_data = trim_all(matrix_to_cell_lines(read_file(file_path)));
table_cells = cellfun(@(x) parse_line(x), table_data, 'UniformOutput', 0);
table_cells = table_cells(2:end, :);

pubmed_table = cell(length(table_cells), 5);

for ii = 1:length(table_cells)
    pubmed_table(ii, :) = table_cells{ii};
end
end

function parsed_line = parse_line(x)
% Parses a line by breaking it into cells based on tab characters.
%
% Input:
%   - x: Input line to be parsed.
%
% Output:
%   - parsed_line: A cell array containing substrings obtained by splitting the input line at tab positions.

tabIndices = find(x == 9);

% If no tab character is found, create indices for the whole line
if isempty(tabIndices)
    tabIndices = [1, length(x) + 1];
end

parsed_line = {x(1:tabIndices(1) - 1), ...
              x(1 + tabIndices(1):tabIndices(2) - 1), ...
              x(1 + tabIndices(2):tabIndices(3) - 1), ...
              x(1 + tabIndices(3):tabIndices(4) - 1), ...
              x(1 + tabIndices(4):end)};
end

function [string_matrix, line_size] = read_file(file_path)
% Reads a text file and returns a matrix of strings.
%
% Inputs:
%   - file_path: The path to the text file to be read.
%
% Outputs:
%   - string_matrix: A matrix of strings where each row corresponds to a line in the file.
%   - line_size: The fixed size of each line in the file.

fid = fopen(file_path, 'r');
file_content = fread(fid, '*char')';
fclose(fid);

% Remove carriage return characters
file_content(file_content == 13) = [];

line_size = max(abs(diff(find(file_content - 0 == 10))));
num_lines = length(find(file_content - 0 == 10));
string_matrix = 32 * ones(num_lines, line_size);

line_breaks = find(file_content - 0 == 10);
line_breaks = [0 line_breaks]; % Corrected on 01/08/2013

for i = 1:(length(line_breaks) - 1)
    string_matrix(i, 1:(line_breaks(i + 1) - line_breaks(i) - 1)) = ...
        file_content((line_breaks(i) + 1):(line_breaks(i + 1) - 1));
end
string_matrix = char(string_matrix);
end

function cell_lines = matrix_to_cell_lines(matrix)
% Converts a matrix into a cell array where each row becomes a cell.
%
% Inputs:
%   - matrix: The input matrix to be converted.
%
% Output:
%   - cell_lines: A cell array with rows of the matrix as individual cells.

cell_lines = mat2cell(matrix, ones(1, size(matrix, 1)), size(matrix, 2));
end

function trimmed_strings = trim_all(xstr, varargin)
% Trim All
%
% This function removes specified characters or spaces from the beginning and
% end of each string in a cell array or a single string.
%
% Parameters:
%   xstr: Input string or cell array of strings to be trimmed
%   varargin: Optional parameter specifying the characters to be trimmed
%
% Output:
%   trimmed_strings: Trimmed string or cell array of trimmed strings

% Set the default trimming character to space (32)
trim_char = ' ';

% Check if varargin is not empty, if yes, use the specified character
if ~isempty(varargin)
    trim_char = varargin{1};
end

% Check if xstr is a cell array
if iscell(xstr)
    % Apply trim_edges to each element in the cell array
    trimmed_strings = cellfun(@(x) trim_edges(x, trim_char), xstr, 'UniformOutput', false);
else
    % Apply trim_edges to the input string
    trimmed_strings = trim_edges(xstr, trim_char);
end
end

function trimmed_str = trim_edges(xstr, ch)
% Trim Leading and Trailing Characters
%
% This function trims leading and trailing occurrences of a specified character
% from the input string.
%
% Parameters:
%   xstr: Input string to be trimmed
%   ch: Character to be trimmed from the edges
%
% Output:
%   trimmed_str: Resulting string after trimming

% Find indices of non-ch characters in the string
non_ch_indices = find(~(xstr == ch));

% Determine the substring between the first and last non-ch characters
trimmed_str = xstr(min(non_ch_indices):max(non_ch_indices));
end
