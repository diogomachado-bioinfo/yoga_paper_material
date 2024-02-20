function [word_tree] = create_global_tree(PUBMED_STRUCT, flag_word, ...
    varargin)
% CREATE_GLOBAL_TREE - Create a global dendrogram of words (word tree) for
% relevant words based on PubMed data.
%
% Input:
%   - PUBMED_STRUCT: PubMed data structure.
%   - flag_word: Target word for tree construction.
%   - varargin: Optional parameter for saving the phylogenetic tree.
%
% Output:
%   - word_tree: Phylogenetic tree structure.
%
% This function selects relevant words based on frequency, relatedness, and
% time correlation, then builds a word tree using the selected
% words and their embeddings.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_word = upper(flag_word);

% Step 1: Filtering words based on frequency and length
process_num = 0;
process_num_max = 2;

% Frequency of words in the order of PUBMED_STRUCT.word_list
yfr0 = PUBMED_STRUCT.normalized_freq(PUBMED_STRUCT.unique_hashes);
inonum = cellfun(@(x) strcmp(x, only_letters(x,1)), ...
    PUBMED_STRUCT.word_list, 'Un', 1);
PUBMED_STRUCT.filter_indices = (yfr0 > 1.5) & (cellfun(@length, ...
    PUBMED_STRUCT.hash_indices) > 1) & (inonum);

% Selecting relevant words based on relatedness
n = length(PUBMED_STRUCT.word_list_filter);
ktop = 20;
iTops = zeros(n, ktop);

process_num = process_num + 1;
reverse_str = '';
fprintf('\t[%02i/%02i] Select relevant words: ', process_num, ...
    process_num_max);
for ii = 1:n
    iiWRD = PUBMED_STRUCT.word_list_filter(ii);
    [~, iidstop] = related_words(iiWRD, PUBMED_STRUCT, 2, ktop);
    iTops(ii,:) = iidstop;

    msg = sprintf('%i/%i', ii, n);
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n');

% Counting words by text
process_num = process_num + 1;
reverse_str = '';
fprintf('\t[%02i/%02i] Count words by text: ', process_num, ...
    process_num_max);
fnucleo = zeros(n, 1);
for ii = 1:n
    xword = PUBMED_STRUCT.word_list_filter(ii);
    [~, u] = search_word(PUBMED_STRUCT, xword, 'k_first', 50);
    fnucleo(ii) = u;
    
    msg = sprintf('%i/%i', ii, n);
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n');

% Calculating scores and preparing for tree construction
crtop = sort_indices(PUBMED_STRUCT.word_time_corr(:,1)); % Consistent
                                                         % growth over time
scocr = 1 - crtop / n;
sconuc = fnucleo / 50;

score3 = 0.05 * scocr + 0.95 * sconuc; % Low weight for count
xscore = score3;
iscore = sort_indices(xscore, 'descend');

WRDcrUp = lower(PUBMED_STRUCT.word_list_filter);
icut= round(0.30*n);
WRDcrUp(iscore(1:icut)) = upper(WRDcrUp(iscore(1:icut)));
WRDcrUp = cell_to_cell_lines([num_to_str_cell(1:n) WRDcrUp ...
    num_to_str_cell(xscore)]);

ntree = min(1500, length(WRDcrUp(:,1))); 

idsn = find(str_find_logical(WRDcrUp, flag_word));
iroot = WRDcrUp{idsn(cellfun(@(x) find_approx_substrings(flag_word, x), ...
    PUBMED_STRUCT.word_list_filter(idsn), 'Un',1))};

% Building the phylogenetic tree
word_tree = generate_phylogenetic_tree(...
    PUBMED_STRUCT.word_embedding_filter_pc(iscore(1:ntree),1:300), ...
    WRDcrUp(iscore(1:ntree)), 'Method', 'Nj', 'Distance', 0.3, ...
    'Reroot', iroot, 'Show', 0);

if numel(varargin)
    phytreewrite(varargin{1}, word_tree);
end

end

function [logical_indices] = str_find_logical(xCell, xexpr)
% STR_FIND_LOGICAL - Search cells that satisfy the logical expression xexpr.
%
% This function performs a logical search on a cell array xCell based on 
% a logical expression xexpr. The logical expression can contain 
% operators such as & (AND), | (OR), ~ (NOT), and parentheses, as well as 
% search strings. Spaces are considered delimiters for operators, and '_' 
% can be used to include spaces in the search.
%
% Input:
%   - xCell: Cell array to be searched.
%   - xexpr: Logical expression for the search.
%
% Output:
%   - logical_indices: Logical values indicating elements in xCell that
%     satisfy xexpr.
%
% Example:
%   logical_indices = str_find_logical(xCell, 'word1 & (word2 | word3)');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preprocess the logical expression
xexpr = strrep(xexpr, '(', ' ( ');
xexpr = strrep(xexpr, ')', ' ) ');
xexpr = strrep(xexpr, '~', ' ~ ');
xexpr = strrep(xexpr, '|', ' | ');
xexpr = strrep(xexpr, '&', ' & ');
xexpr = eliminate_repetitions(xexpr);

% Convert the logical expression into a cell array of words
p2w = phrase_to_words({xexpr});
ichan = cellfun(@(x) max(ismember(x, '&|~()')), p2w, 'Un', 1);
p2w(ichan) = cellfun(@(x) [' ' trim_all(x) ' '], p2w(ichan), 'Un', 0);
p2w(~ichan) = cellfun(@(x) ['str_find_cell(xCell, #' x '#)'], ...
    p2w(~ichan), 'Un', 0);

% Eliminate repetitions and convert to a string
xstr = eliminate_repetitions(char(mat_to_vec(char(p2w))));
xstr(xstr=='#') = '''';

% Evaluate the logical expression to obtain logical values
logical_indices = eval(xstr)';  % Logical values - Elements pointed to
                                % in xCell
end

function [similarity_scores, matching_indices] = find_approx_substrings(...
    x_str, substring, varargin)
% Find approximately matching substrings of 'substring' in 'x_str'
%
% x_str - Input string
% substring - Substring to find in 'x_str'
% varargin - Optional parameter for the similarity cutoff

substring_length = length(substring);

% Convert the input string to a list of substrings of the same length as
% the target substring
substring_list = dna_to_list(x_str, substring_length);
list_length = length(substring_list(:, 1));

% Calculate similarity scores by comparing each substring in the list with
% the target substring
similarity_scores = sum((substring_list == repmat(substring, ...
    list_length, 1)), 2) / substring_length;

% Check if a similarity cutoff is provided in varargin
if isempty(varargin)
    matching_indices = find(similarity_scores == 1);
else
    similarity_cutoff = varargin{1};
    matching_indices = find(similarity_scores >= similarity_cutoff);
end

% Extract the similarity scores for the matching indices
similarity_scores = similarity_scores(matching_indices);

% If no matches are found, set similarity_score to 0
if isempty(matching_indices)
    similarity_scores = 0;
end
end

function patterns_matrix = dna_to_list(dna_sequence, pattern_length)
% Create a matrix of patterns (subsequences) for a given DNA sequence
%
% dna_sequence - Input DNA sequence
% pattern_length - Length of the patterns to extract

sequence_length = length(dna_sequence);
column_indices = uint32(1:(sequence_length - pattern_length + 1))';

% Create matrix indices for extracting patterns of the specified length
pattern_indices = repmat(column_indices, 1, pattern_length) + ...
    repmat(uint32(1:pattern_length), (sequence_length - ...
    pattern_length + 1), 1) - 1;

% Extract patterns from the DNA sequence using the calculated indices
patterns_matrix = dna_sequence(pattern_indices);

% If the pattern length is 1, transpose the matrix to have patterns in rows
if pattern_length == 1
    patterns_matrix = patterns_matrix';
end
end