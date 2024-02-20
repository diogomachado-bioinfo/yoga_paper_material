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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the size of the input vector
vector_size = size(input_vector);

% Initialize the output as the input vector
column_vector = input_vector;

% Check if the vector is a row vector, then transpose it to make it a column vector
if vector_size(2) > vector_size(1)
    column_vector = input_vector';
end

end
