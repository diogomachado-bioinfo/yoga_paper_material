function [mret, idstop, otree] = related_words(input_words, ...
    PUBMED_STRUCT, cut_min, k_list, varargin)
% Searches for related words in the filtered list.
%
% Input:
%   - input_words: Words to search for, can be a single word or a cell array.
%   - PUBMED_STRUCT
%   - cut_min: Minimum occurrences of words for inclusion.
%   - k_list: Number of relevant words to return.
%   - varargin: Additional parameters (optional).
%
% Output:
%   - mret: Matrix of relevant words.
%   - idstop: Indices of the stop words.
%   - otree: Word tree based on relevant words.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = varargin;
m_show = find_in_varargin(V, 'Show', 0);
n_tree = find_in_varargin(V, 'ntree', 30);
X = upper(input_words);

Wwords = PUBMED_STRUCT.word_embedding_filter_pc(:,1:50);
Wwords0 = PUBMED_STRUCT.word_embedding_pc(:,1:50); % Também importante

X_WORD = find(ismember(PUBMED_STRUCT.unique_hashes, words_to_hashinrecs(X)));
ids = sort_indices(norm_vect(repmat(mean(Wwords0(X_WORD,:), 1), ...
    length(PUBMED_STRUCT.word_list_filter), 1) - Wwords), 'descend');

if m_show
    wlist = indexa(PUBMED_STRUCT.word_list_filter(ids), (cellfun(@length, PUBMED_STRUCT.word_list_filter(ids), ...
        'Un', 1) > 1) & cellfun(@(x) length(x) > cut_min, PUBMED_STRUCT.hash_indices_filter(ids), ...
        'Un', 1));
    disp(['Search base: ' PUBMED_STRUCT.word_list(X_WORD)'])
else
    wlist = index_elements(PUBMED_STRUCT.word_list_filter(ids), (cellfun(@length, ...
        PUBMED_STRUCT.word_list_filter(ids), 'Un', 1) > 1) & cellfun(@(x) length(x) > ...
        cut_min, PUBMED_STRUCT.hash_indices_filter(ids), 'Un', 1));
    PUBMED_STRUCT.word_list(X_WORD);
end

wbase = PUBMED_STRUCT.word_list(X_WORD);
end_x = length(wlist);
idsshow = 1 - k_list + end_x:end_x;
if m_show
    mret = cell2mat(mat_to_vec([(mat_to_vec(wlist(idsshow))) ...
        mat_to_cell_lines(create_space_matrix(k_list, 1))]));
    disp(mret);
else
    mret = cell2mat(mat_to_vec([(mat_to_vec(wlist(idsshow))) ...
        mat_to_cell_lines(create_space_matrix(k_list, 1))]));
end
mret = eliminate_repetitions(strrep(mret, PUBMED_STRUCT.word_list{X_WORD(1)}, ''));
idstop = index_elements(flipud(ids), 1:k_list);
xx = index_elements(flipud(ids), 1:n_tree);

if m_show
    otree = generate_phylogenetic_tree(Wwords(xx, :), PUBMED_STRUCT.word_list_filter(xx), ...
        'Method', 'Nj', 'Distance', 0.3, 'Reroot', wbase);
else
    otree = generate_phylogenetic_tree(Wwords(xx, :), PUBMED_STRUCT.word_list_filter(xx), ...
        'Method', 'Nj', 'Distance', 0.3, 'Reroot', wbase, 'Show', 0);
end
end

function selected_elements = index_elements(X, indices_or_conditions)
% Selects elements from X based on specified indices or logical conditions.
%
% Input:
%   - X: Input vector or array.
%   - indices_or_conditions: Indices or logical conditions for selection.
%
% Output:
%   - selected_elements: Selected elements from X.

if ischar(indices_or_conditions)
    % If indices_or_conditions is a string, evaluate it as a logical condition.
    selected_elements = eval(['X(' indices_or_conditions ')']);
else
    % If indices_or_conditions is numeric, directly index the elements.
    selected_elements = X(indices_or_conditions);
end
end

function space_matrix = create_space_matrix(rows, columns)
% Creates a matrix filled with space characters (ASCII 32).
%
% Input:
%   - rows: Number of rows in the matrix.
%   - columns: Number of columns in the matrix.
%
% Output:
%   - space_matrix: Matrix filled with space characters.

space_matrix = char(32 * ones(rows, columns));
end
