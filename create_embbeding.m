function PUBMED_STRUCT = create_embbeding(PUBMED_STRUCT)
% CREATE_EMBBEDING - Create word and text embeddings for PubMed data.
%
% This function takes a structure PUBMED_STRUCT containing PubMed data 
% and calculates word and text embeddings based on the provided data. 
% The resulting embeddings are stored back in the PUBMED_STRUCT structure.
%
% Input:
%   - PUBMED_STRUCT: Structure containing PubMed data.
%
% Output:
%   - PUBMED_STRUCT: Updated structure with word and text embeddings.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize process variables
process_num = 0;
process_num_max = 2;

% Part 1: Creating Word Embedding
word_embedding = zeros(PUBMED_STRUCT.num_unique_words, ...
    PUBMED_STRUCT.num_sweep_rows);
process_num = process_num + 1;
reverse_str = '';
fprintf('\t[%02i/%02i] Creating word embedding: ', process_num, ...
    process_num_max);
for ii = 1:PUBMED_STRUCT.num_unique_words
    iinds = PUBMED_STRUCT.hash_indices{ii};
    word_embedding(ii,:) = mean(PUBMED_STRUCT.sweep_result(iinds,:));
    
    % Display progress
    msg = sprintf('%i/%i', ii, PUBMED_STRUCT.num_unique_words);
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n');

% Apply principal component analysis (PCA) on the word embedding
word_embedding_filter = word_embedding(PUBMED_STRUCT.filter_indices,:);
word_embedding_filter_pc = word_embedding_filter * princomp(...
    word_embedding_filter);
% To have the entire list of words for searching
word_embedding_pc = word_embedding * princomp(word_embedding_filter);

% Part 2: Creating Text Embedding
text_embedding = zeros(PUBMED_STRUCT.num_texts, length(...
    word_embedding(1,:)));
iTEXTcln = cell(PUBMED_STRUCT.num_texts, 1);
process_num = process_num + 1;
reverse_str = '';
fprintf('\t[%02i/%02i] Creating text embedding: ', process_num, ...
    process_num_max);
for ii = 1:PUBMED_STRUCT.num_texts
    iwrd = PUBMED_STRUCT.index_table(PUBMED_STRUCT.text_hashes{ii});
    iTEXTcln{ii} = iwrd;
    text_embedding(ii,:) = mean(word_embedding(iwrd,:), 1);
    
    % Display progress
    msg = sprintf('%i/%i', ii, PUBMED_STRUCT.num_texts);
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n');

% Update the structure with new embeddings
PUBMED_STRUCT.word_embedding = word_embedding;
PUBMED_STRUCT.word_embedding_filter = word_embedding_filter;
PUBMED_STRUCT.word_embedding_filter_pc = word_embedding_filter_pc;
PUBMED_STRUCT.word_embedding_pc = word_embedding_pc;
PUBMED_STRUCT.text_embedding = text_embedding;

end
