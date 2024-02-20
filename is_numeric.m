function result = is_numeric(C)
% Checks if elements in the input are numeric.
%
% Input:
%   - C: Input cell array or element.
%
% Output:
%   - result: True for numeric elements.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    result = ~ischar(C) & ~iscell(C);

    if iscell(C)
        result = cellfun(@is_numeric, C, 'UniformOutput', 1);
    end
end
