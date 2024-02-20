function [search_results, word_occurrences] = search_word(PUBMED_STRUCT, ...
    x_word, varargin)
% SEARCH_WORD - Search for a word in articles from PUBMED_STRUCT.
%
% This function searches for a specific word in the text of PubMed articles 
% stored in the provided PUBMED_STRUCT. It then presents the k_first 
% top-ranked articles that contain the searched word. The function also 
% returns the number of occurrences of the word in the presented articles.
%
% Input:
%   - x_word: Word to search for.
%   - PUBMED_STRUCT: Structure containing PubMed data.
%   - varargin: Optional parameter-value pairs for customization.
%       - 'custom_text': Array of texts corresponding to the processed
%         articles (default: PUBMED_STRUCT.processed_text). NOTE: The
%         'custom_text' must be in correspondence with
%         PUBMED_STRUCT.processed_text.
%       - 'k_first': Number of articles presented in the search
%         (default: 20).
%
% Output:
%   - search_results: Relevant information from the search.
%   - word_occurrences: Number of occurrences of the word in the presented
%     articles.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = varargin;
custom_text = find_in_varargin(V, 'custom_text', ...
    PUBMED_STRUCT.processed_text);
k_first = find_in_varargin(V, 'k_first', 20);

% Ensure k_first is within bounds
k_first = min(length(custom_text), k_first);

% Number of texts and word indices for the searched word
num_texts = length(PUBMED_STRUCT.text_embedding(:, 1));
word_indices = PUBMED_STRUCT.index_table(words_to_hashinrecs(upper(...
    x_word)));

% Extract unique words and calculate mean word embedding for the searched
% word
PUBMED_STRUCT.word_list(word_indices); 
word_embedding_mean = mean(PUBMED_STRUCT.word_embedding(...
    word_indices, :), 1);

% Rank articles based on similarity to the mean word embedding
ranked_indices = sort_indices(norm_vect(...
    PUBMED_STRUCT.text_embedding(:, :) - repmat(...
    word_embedding_mean, num_texts, 1)));

% Take the first index for analysis
chosen_index = ranked_indices(1);

% Extract information about terms present in the top-ranked articles
word_info = cellfun(@(x) PUBMED_STRUCT.index_table(unique(...
    PUBMED_STRUCT.text_hashes{chosen_index}(ismember(...
    PUBMED_STRUCT.text_hashes{chosen_index}, x)))), ...
    PUBMED_STRUCT.text_hashes(ranked_indices(:)), 'Un', 0);

% Count occurrences of terms in the top-ranked articles
word_counts = cellfun(@(x) length(x), word_info);

% Create a cell array with information about the top-ranked articles
search_results = cell_to_cell_lines([num_to_str_cell(1:k_first) ...
    num_to_str_cell(ranked_indices(1:k_first)) num_to_str_cell(...
    word_counts(1:k_first)) custom_text(ranked_indices(1:k_first))]);

% Count occurrences of the searched word in the presented articles
occurrences_with_word = str_find_cell(search_results, x_word{1});
word_occurrences = sum(occurrences_with_word);

end
