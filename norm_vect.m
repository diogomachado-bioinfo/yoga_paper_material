function mret = norm_vect(mat)
% Calculates the L2 norm for each row of the matrix.
%
% Input:
%   - mat: Input matrix.
%
% Output:
%   - mret: Vector containing L2 norms for each row.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mret = sqrt(sum(mat.^2, 2));
end
