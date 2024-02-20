function column_vector = vec_in_col(vec)
% Ensures the input vector is always a column vector.
%
% Input:
%   - vet: Input vector.
%
% Output:
%   - column_vector: Column vector.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = size(vec);
column_vector = vec;

if (s(2) > s(1))
    column_vector = vec';
end    
end
