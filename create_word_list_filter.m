function PUBMED_STRUCT = create_word_list_filter(PUBMED_STRUCT)
% CREATE_WORD_LIST_FILTER - Filters words in PUBMED_STRUCT based on
% specified criteria.
%
% This function filters words in the given PUBMED_STRUCT based on certain
% criteria, including word frequencies, presence of numbers, and additional
% criteria related to hash collisions. It updates the structure with the
% filtered word list and associated information.
%
% Input:
%   - PUBMED_STRUCT: Structure containing PubMed data.
%
% Output:
%   - PUBMED_STRUCT: Updated structure with filtered word list and related
%     information.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate word frequencies
word_freqs = PUBMED_STRUCT.normalized_freq(PUBMED_STRUCT.unique_hashes);

% Identify indices of words containing only letters
no_num_idx = cellfun(@(x) strcmp(x, only_letters(x, 1)), ...
    PUBMED_STRUCT.word_list, 'Un', 1);
filter_indices = (word_freqs > 1.5) & ...
    (cellfun(@length, PUBMED_STRUCT.hash_indices) > 1) & (no_num_idx);
filter_indices_find = find(filter_indices);
word_list_filter = PUBMED_STRUCT.word_list(filter_indices);

% Resolve hash collisions
nwrd = length(word_list_filter);
idemword = zeros(length(filter_indices), 1);
yfr = word_freqs(filter_indices);
hash_indices_filter = PUBMED_STRUCT.hash_indices(filter_indices);
for ii = 1:nwrd
    % Check if the word is in the dictionary and meets additional criteria
    isindic = strcmp(PUBMED_STRUCT.RSMFRQ.WRDIRdic(...
        words_to_hashinrecs(word_list_filter{ii})), word_list_filter{ii});
    if ~isindic && (yfr(ii) > 5) && (length(hash_indices_filter{ii}) > 20)
        isindic = 1;
    end
    idemword(filter_indices_find(ii)) = isindic;
end

% Create variables with the filter
filter_indices(~logical(idemword)) = false;
word_list_filter = PUBMED_STRUCT.word_list(filter_indices);
hash_indices_filter = PUBMED_STRUCT.hash_indices(filter_indices);
unique_hashes_filter = PUBMED_STRUCT.unique_hashes(filter_indices);

% Synchronize HASHES with PUBMED_STRUCT.filter_indices
xmod = length(PUBMED_STRUCT.word_hash);
index_table_filter = zeros(xmod, 1);
index_table_filter(unique_hashes_filter) = 1:length(word_list_filter);

% Insert new data into the struct
PUBMED_STRUCT.word_list_filter = word_list_filter;
PUBMED_STRUCT.hash_indices_filter = hash_indices_filter;
PUBMED_STRUCT.unique_hashes_filter = unique_hashes_filter;
PUBMED_STRUCT.filter_indices = filter_indices;
PUBMED_STRUCT.index_table_filter = index_table_filter;

end
