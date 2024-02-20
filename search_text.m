function [search_results] = search_text(PUBMED_STRUCT, ...
    search_substring_or_index, varargin)
% SEARCH_TEXT - Search for a text substring in a PUBMED_STRUCT.
%
% This function uses text embeddings to find articles similar to a given
% text substring or indexed paper in the PUBMED_STRUCT. It calculates the
% similarity between the searched text and all other articles, then
% presents the k_first top-ranked articles based on this similarity
% measure. By default, the output is constructed using 
% PUBMED_STRUCT.processed_text, but it can be replaced by the custom_text
% optional input using the corresponding parameter.
%
% INPUTS:
%   - PUBMED_STRUCT: Structure containing processed information from PUBMED
%     articles.
%   - search_substring_or_index: Substring of the text to search for in the
%     processed articles or index of the paper.
%   - varargin: Optional parameter-value pairs for customization.
%       - 'custom_text': Array of texts corresponding to the processed
%         articles (default: PUBMED_STRUCT.processed_text). NOTE: The
%         'custom_text' must be in correspondence with
%         PUBMED_STRUCT.processed_text.
%       - 'k_first': Number of top-ranked articles to present in the search
%         results (default: 20).
%
% OUTPUT:
%   - search_results: Resulting cell array containing information about the
%     top-ranked articles.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = varargin;
custom_text = find_in_varargin(V, 'custom_text', ...
    PUBMED_STRUCT.processed_text);
k_first = find_in_varargin(V, 'k_first', 20);

% Ensure k_first is within bounds
k_first = min(length(custom_text), k_first);

% Number of texts and check if searching by substring or index
num_texts = length(PUBMED_STRUCT.text_embedding(:, 1));
if ischar(search_substring_or_index) % Search by substring
    search_pattern = search_substring_or_index;
    indices_to_take = strfindcell(PUBMED_STRUCT.processed_text, ...
        upper(search_pattern));
else % Search by index
    indices_to_take = search_substring_or_index;
end

% Take the first index if multiple indices are found
index_to_take = indices_to_take(1);

% Extract the embedding of the searched text
search_embedding = PUBMED_STRUCT.text_embedding(index_to_take,:);

% Calculate similarity and rank articles based on similarity
ranked_indices = sort_indices(norm_vect(...
    PUBMED_STRUCT.text_embedding(:,:) - ...
    repmat(search_embedding, num_texts, 1)));

% Extract information about terms present in the top-ranked articles
top_words_info = cellfun(@(x) PUBMED_STRUCT.index_table(unique(...
    PUBMED_STRUCT.text_hashes{index_to_take}(ismember(...
    PUBMED_STRUCT.text_hashes{index_to_take}, x)))), ...
    PUBMED_STRUCT.text_hashes(ranked_indices(:)), 'Un', 0);
word_counts = cellfun(@(x) length(x), top_words_info);

% Create a cell array with information about the top-ranked articles
search_results = cell_to_cell_lines([num_to_str_cell(1:k_first) ...
    num_to_str_cell(ranked_indices(1:k_first)) num_to_str_cell(...
    word_counts(1:k_first)) custom_text(ranked_indices(1:k_first))]);

end
