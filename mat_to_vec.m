function vector = mat_to_vec(matrix)
% Converts the matrix into a vector, by rows.
%
% Input:
%   - matrix_input: Input matrix.
%
% Output:
%   - vector: Vector obtained by rearranging the elements of the matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = matrix';
n = size(q);
if ~iscell(matrix)
    x = zeros(1, prod(n));
    x(1:prod(n)) = q(q | ~q)';
else
    qmat = ones(n(1), n(2));
    x = q(qmat == 1)';
end
vector = x;
end
