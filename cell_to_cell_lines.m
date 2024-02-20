function cell_lines = cell_to_cell_lines(C)
% Converts a matrix of cells with multiple columns into cell lines.
%
% Input:
%   - C: Matrix of cells with multiple columns.
%
% Output:
%   - cellLines: Cell array containing lines formed by concatenating 
%     elements of each row in the matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_lines = {'Empty'};
if ~isempty(C)
    [~, m] = size(C);
    C_ = mat_to_vec(C)';
    idsChange = find(is_numeric(C_));
    
    if ~isempty(idsChange)
        C_{is_numeric(C_)} = cellfun(@num2str, C_(idsChange), 'Un', 1);
    end
    
    C_ = vec_in_col(cellfun(@(x) [x ' '], C_, 'Un', 0)');
    C2 = vec_to_mat(C_, m);
    
    [rows, ~] = size(C2);
    if rows > 1
        cell_lines = trim_all(cellfun(@(x) [C2{x, :}], mat_to_cell_lines((1:rows)'), 'Un', 0));
    else
        cell_lines = {eliminate_repetitions(char(mat_to_vec(char(cellfun(@(x) [x ' '], C2, 'Un', 0)))))};
    end
end
end

function matrix = vec_to_mat(vector, width)
% Creates a matrix with the values of the input vector by row with the specified width.
% If the last row is incomplete, it is filled with zeros.
%
% Input:
%   - vector: Input vector.
%   - width: Width of the matrix columns.
%
% Output:
%   - matrix: Output matrix.

num_values = length(vector);
num_rows = ceil(num_values / width);
indices = find(ones(width, num_rows));

if iscell(vector)
    matrix = cell(width, num_rows);
else
    matrix = zeros(width, num_rows);
end

matrix(indices(1:num_values)) = vector;
matrix = matrix';

end
