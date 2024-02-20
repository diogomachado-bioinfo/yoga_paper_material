function modified_sequence = eliminate_repetitions(seq, varargin)
% Eliminate Repetitions
%
% This function removes consecutive occurrences of a specified character or
% the space character (32) from the input sequence.
%
% Parameters:
%   seq: Input sequence
%   varargin: Optional parameter specifying the character to be removed
%
% Output:
%   modified_sequence: Sequence with only one occurrence of the specified character
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if varargin is empty, if not, use the specified character; otherwise, use space (32)
if isempty(varargin)
    character = ' ';
else
    character = varargin{1};
end

% Check if the input sequence is a cell array
if iscell(seq)
    % Apply make_one_occurrence to each element of the cell array
    modified_sequence = cellfun(@(x) make_one_occurrence(x, character), seq, 'UniformOutput', false);
else
    % Apply make_one_occurrence to the input sequence
    modified_sequence = make_one_occurrence(seq, character);
end

end

function modified_sequence = make_one_occurrence(seq, character)
% Make One Occurrence
%
% This function removes consecutive occurrences of a specified character in a sequence.
%
% Parameters:
%   seq: Input sequence
%   character: Character to be removed from consecutive occurrences
%
% Output:
%   modified_sequence: Sequence with only one occurrence of the specified character

% Create a copy of the input sequence
modified_sequence = seq;

% Remove consecutive occurrences of the specified character
occurrence_indices = strfind(modified_sequence, [character character]);
modified_sequence(occurrence_indices) = [];

% Trim all occurrences of the specified character
modified_sequence = trim_all(modified_sequence, character);

end
