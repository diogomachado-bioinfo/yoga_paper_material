function PUBMED_STRUCT = temporal_correlation(PUBMED_STRUCT)
% TEMPORAL_CORRELATION - Analyzes the temporal correlation of word
% occurrences in PubMed data.
%
% This function calculates the temporal correlation of word occurrences in
% PubMed articles over time. It provides insights into how the frequency
% of specific words correlates with the passage of time. The analysis
% includes calculating correlation coefficients, p-values, and storing
% results for further exploration.
%
% Input:
%   - PUBMED_STRUCT: Structure containing PubMed data.
%
% Output:
%   - PUBMED_STRUCT: Updated structure with additional temporal correlation
%     information.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the starting year for analysis as one year before the minimum
% publication year
anocut = min(PUBMED_STRUCT.publication_years) - 1;

% Set the ending year for analysis as the maximum publication year
anofim = max(PUBMED_STRUCT.publication_years);   % Generalize

% Count occurrences of publication years after the cut-off year
ctano = count_occurrences(PUBMED_STRUCT.publication_years(...
    PUBMED_STRUCT.publication_years > anocut));

% Calculate the length of the filtered word list in PUBMED_STRUCT
nwrd = length(PUBMED_STRUCT.word_list_filter);

% Initialize an array to store a matrix to correlation statistics for each
% word, depicting the relationship between the occurrence of articles per
% year and the passage of time. Each row corresponds to a specific word,
% with the first column representing the correlation coefficient, and the
% second column containing the associated p-value. These statistics provide
% insights into the temporal correlation patterns associated with each word
% in the dataset.
word_time_corr = zeros(nwrd, 2);

% Calculate the number of years for analysis
nanos = anofim - anocut;

% Initialize an array to store the list of years and their counts
count_by_year = zeros(nanos, 1);

% Update the list of years with counts from the calculated occurrences
iano = (ctano(:, 1) - anocut); 
iano = iano(iano > 0);
count_by_year(iano) = ctano(:, 2);
count_by_year = [(1 + anocut:anofim)' count_by_year];

% Initialize an array to store a matrix representing word occurrences in
% texts for each year. Each row corresponds to a specific word, and each
% column represents the count of occurrences of that word in the respective
% year. The matrix provides a detailed record of word occurrences over
% time, allowing analysis of temporal patterns in the corpus.
word_occ_by_year = zeros(nwrd, nanos);

% Loop through each word in the filtered word list
for ii = 1:nwrd
    % Initialize an array to store the list of occurrences for each year
    lista = zeros(nanos, 1);

    % Find indices of occurrences of the current word in processed text
    idsano = str_find_cell(PUBMED_STRUCT.processed_text, ...
        PUBMED_STRUCT.word_list_filter{ii})';

    % Check if there are any occurrences
    if sum(idsano) > 0 
        % Count occurrences of publication years after the cut-off year
        % for the current word
        ctii = count_occurrences(PUBMED_STRUCT.publication_years(...
            idsano & PUBMED_STRUCT.publication_years > anocut));

        % Update the list of occurrences with counts from the calculated
        % occurrences
        iano = (ctii(:, 1) - anocut);
        lista(iano) = ctii(:, 2);

        % Calculate the correlation with the overall counts
        xcr = lista ./ max(count_by_year(:, 2), 1);

        % Apply fuzzy moving average to smooth the data
        xcr = fuzzy_moving_average(xcr, 1, [0.9, 2])';

        % Store the results in the word_occ_by_year array
        word_occ_by_year(ii, :) = lista';

        % Calculate correlation and p-value for the smoothed data
        xycorr = [xcr(1:end-1), (1:nanos-1)'];
        if ~isempty(xycorr)
            [cr, p] = corr(xycorr(:, 1), xycorr(:, 2));
        else
            cr = 0;
            p = 1;
        end
    else % Word does not occur in the period
        cr = 0;
        p = 1;
    end

    % Store correlation and p-value in the word_time_corr array
    word_time_corr(ii, :) = [cr p];
    
    PUBMED_STRUCT.count_by_year = count_by_year;
    PUBMED_STRUCT.word_occ_by_year = word_occ_by_year;
    PUBMED_STRUCT.word_time_corr = word_time_corr;
end
end