function words_cell = phrase_to_words(cell_of_strings)
% Phrase to Words
%
% This function takes a cell array of strings and returns a cell array with
% each word separated in the input strings.
%
% Parameters:
%   cell_of_strings: Cell array of strings to be processed
%
% Output:
%   words_cell: Cell array with individual words extracted from the input strings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply split_string to each element in the cell array
split_strings_cell = cellfun(@(x) split_string(x), cell_of_strings, 'UniformOutput', false);

% Combine the resulting cell array into a single cell array of words
words_cell = [split_strings_cell{:}];

end
