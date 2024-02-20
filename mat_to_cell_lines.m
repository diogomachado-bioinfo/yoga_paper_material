function cell_lines = mat_to_cell_lines(matrix)
% Generates a cell array with lines from the input matrix.
%
% Input:
%   - matrix: Input matrix.
%
% Output:
%   - cell_lines: Cell array containing lines of the matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_lines = mat2cell(matrix, ones(1, length(matrix(:,1))),length(matrix(1,:)));
end
