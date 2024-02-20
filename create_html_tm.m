function [output_dir] = create_html_tm(PUBMED_STRUCT, flag_word, varargin)
% Create HTML files based on PubMed data.
%
% Parameters:
%   PUBMED_STRUCT: Structure containing PubMed data.
%   flag_word: String used in the output directory name.
%   varargin: Optional parameter-value pairs for customization.
%       - 'words_html': Generate HTML files related to words (default:
%         true).
%       - 'text_html': Generate HTML files related to text data (default:
%         true).
%
% Returns:
%   output_dir: Directory where HTML files are created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract optional variables from varargin
V = varargin;
words_html = find_in_varargin(V, 'words_html', true);
text_html = find_in_varargin(V, 'text_html', true);

% Create output directory based on the provided flag_word and date
output_dir = [flag_word '_' strrep(date, '-', '_')];

% Initialize variables for tracking the progress of processes
process_num = 0;
process_num_max = 0;

% Check if words_html and text_html are true to update process_num_max
if words_html
    process_num_max = process_num_max + 1;
end
if text_html
    process_num_max = process_num_max + 1;
end

% Create output directory if it does not exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Create DATA directory within the output directory if it does not exist
DATA_dir = fullfile(output_dir, 'DATA');
if ~exist(DATA_dir, 'dir')
    mkdir(DATA_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process for generating HTML files related to words
if words_html
    % Initialize HTML table for words
    CHTML = {'<table border="1">'};
    CHTML{2} = '    <tr>';
    CHTML{3} = '        <td>Cod</td>';
    CHTML{4} = '        <td>WORD</td>';
    CHTML{5} = '        <td>Related words</td>';
    CHTML{6} = '        <td>Link Title</td>';
    CHTML{7} = '        <td>Link Abstracts</td>';
    CHTML{8} = '        <td>Link Graphic</td>';
    CHTML{9} = '    </tr>';

    TX0 = cellfun(@(x, y)  [x '. ' y], ...
        PUBMED_STRUCT.processed_table(:,2), ...
        PUBMED_STRUCT.processed_table(:,3), 'Un', 0);
    TX0 = strrep(upper(TX0),char(9),'');
    TX0 = strrep(upper(TX0),char(10),'');
    TXAN = cellfun(@(x,y) [x char(9) y],TX0,cellfun(@(x) num2str(x), ...
        mat_to_cell_lines(PUBMED_STRUCT.publication_years),'Un',0),'Un',0);

    set(0,'DefaultFigureVisible','off');

    % Preparing WORD.html data
    io = 1;
    process_num = process_num + 1;
    reverse_str = '';
    fprintf('\t[%02i/%02i] Preparing WORDS.html data: ', process_num, ...
        process_num_max);
    max_i = length(PUBMED_STRUCT.word_list_filter);
    for ii=io:max_i
        XWORD = PUBMED_STRUCT.word_list_filter(ii);
        [xstr, ~, otree] = related_words(XWORD, PUBMED_STRUCT, 2, 10);
        xcell = search_word(PUBMED_STRUCT, upper(XWORD), 'custom_text', ...
            TXAN);
        T = char(cellfun(@(x) [x char([9 10])],xcell,'Un',0));

        ctword = [PUBMED_STRUCT.count_by_year(:,1) ...
            (PUBMED_STRUCT.word_occ_by_year(...
            PUBMED_STRUCT.index_table_filter(words_to_hashinrecs(...
            upper(XWORD{1}))),:)'./max(...
            PUBMED_STRUCT.count_by_year(:,2),1))];

        tree_fig_file = fullfile(DATA_dir,[XWORD{1} '.png']);
        year_fig_file = fullfile(DATA_dir,[XWORD{1} '.jpg']);
        abstracts_txt_file = fullfile(DATA_dir,[XWORD{1} '.txt']);

        figure;
        plot(otree);
        saveas(gcf, tree_fig_file);
        close;
        figure;
        hold on;

        if length(PUBMED_STRUCT.count_by_year(:,1))>1
            plot(PUBMED_STRUCT.count_by_year(1:end,1), ...
                fuzzy_moving_average(...
                PUBMED_STRUCT.count_by_year(1:end,2)/...
                sum(PUBMED_STRUCT.count_by_year(1:end,2)),3,[0.9,2]));
            plot(ctword(1:end,1),fuzzy_moving_average(ctword(1:end,2)/...
                sum(ctword(1:end,2)),3,[0.9,2]),'m');
        else
            plot(0,0,'.'); % No year data
        end
        saveas(gcf, year_fig_file);
        close;
        %
        write_file(T, abstracts_txt_file);
        % Put line
        cnt = length(CHTML);
        CHTML{cnt + 1} = '<tr>';
        CHTML{cnt + 2} = ['<td>' num2str(ii) '</td>'];
        CHTML{cnt + 3} = ['<td>' XWORD{1} '</td>'];
        CHTML{cnt + 4} = ['    <td>' lower(xstr) '</td>'];
        CHTML{cnt + 5} = ['    <td><a href="file:.\DATA\' XWORD{1} ...
            '.txt">Abstracts</a></td>'];
        CHTML{cnt + 6} = ['    <td><a href="file:.\DATA\' XWORD{1} ...
            '.png">Tree</a></td>'];
        CHTML{cnt + 7} = ['    <td><a href="file:.\DATA\' XWORD{1} ...
            '.jpg">Year Usage</a></td>'];
        CHTML{cnt + 8} = '</tr>';

        msg = sprintf('%i/%i', ii, max_i);
        fprintf([reverse_str, msg]);
        reverse_str = repmat(sprintf('\b'), 1, length(msg));
    end
    fprintf('\n');

    CHTML{end+1} = '</table>';
    html_file = fullfile(output_dir, 'WORDS.html');
    write_file(char(CHTML), html_file)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process for generating HTML files related to text data
if text_html
    % Initialize HTML table for text data
    CHTML = {'<table border="1">'};
    CHTML{2} = '    <tr>';
    CHTML{3} = '        <td>Cod</td>';
    CHTML{4} = '        <td>Title+Abstracts</td>';
    CHTML{5} = '    </tr>';
    
    TXAN = cellfun(@(x,y) [x char(9) y], PUBMED_STRUCT.processed_text, ...
        cellfun(@(x) num2str(x), mat_to_cell_lines(...
        PUBMED_STRUCT.publication_years),'Un',0),'Un',0);

    process_num = process_num + 1;
    reverse_str = '';
    fprintf('\t[%02i/%02i] Preparing TEXTS.html data: ', process_num, ...
        process_num_max);
    max_i = length(PUBMED_STRUCT.processed_text);
    for ii=1:length(PUBMED_STRUCT.processed_text)
        XTEXT = upper(trim_all(PUBMED_STRUCT.processed_table{ii, 2}));
        xcell = search_text(PUBMED_STRUCT, ii, 'custom_text', TXAN);
        T = char(cellfun(@(x) [x char([9 10])],xcell,'Un',0));
        %
        XNAME = ['ABS' str_zero(ii,5)];
        abstracts_txt_file = fullfile(DATA_dir, [XNAME '.txt']);
        write_file(T, abstracts_txt_file)
        % Put line
        cnt = length(CHTML);
        CHTML{cnt + 1} = '<tr>';
        CHTML{cnt + 2} = ['    <td>' num2str(ii) '</td>'];
        CHTML{cnt + 3} = ['    <td><a href="file:.\DATA\' XNAME ...
            '.txt">' XTEXT '</a></td>'];
        CHTML{cnt + 4} = '</tr>';

        msg = sprintf('%i/%i', ii, max_i);
        fprintf([reverse_str, msg]);
        reverse_str = repmat(sprintf('\b'), 1, length(msg));
    end
    fprintf('\n');

    CHTML{end+1} = '</table>';   
    html_file = fullfile(output_dir, 'TEXTS.html');
    write_file(char(CHTML), html_file)
end
end

function padded_string = str_zero(vectors, padding_length)
% Convert numeric vectors to strings, pad with leading zeros, and return
% the result.
%
% Parameters:
%   - vectors: Numeric vectors to be converted to strings and padded with
%     zeros.
%   - padding_length: Desired length of the resulting string.
%
% Returns:
%   - padded_string: String matrix with leading zeros, representing the
%     input numeric vectors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size_vectors = size(vectors);

% Ensure vectors are treated as columns
if size_vectors(2) > size_vectors(1)
    vectors = vectors';
end

% Create a string matrix with leading spaces and concatenate with numeric
% vectors
string_matrix = [char(32 * ones(length(vectors(:, 1)), padding_length)) ...
    num2str(vectors)];

% Replace spaces with zeros
string_matrix(string_matrix == ' ') = '0';

% Extract the substring with the desired length from the end of each string
padded_string = string_matrix(:, max((end - padding_length + 1), 1):end);

end
