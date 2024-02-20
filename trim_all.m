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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
