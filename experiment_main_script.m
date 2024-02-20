%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WHAT ARE WE LEARNING WITH YOGA: A TEXT MINING APPROACH TO LITERATURE - 
% EXPERIMENT SCRIPT
% -------------------------------------------------------------------------
% Description: Script for run the experiment described in the article "What
% are we learning with Yoga: a text mining approach to literature".
% -------------------------------------------------------------------------
% Authors:
% - Roberto T. Raittz (Code Originality)
% - Diogo de J. S. Machado (Code Organization)
% -------------------------------------------------------------------------
% Date: 02/19/2024
% -------------------------------------------------------------------------
% Note: To run this script, it is necessary to download a required file
%       manually at: <https://bit.ly/4bJRazk>. The script expects this
%       file to be available at current folder or MATLAB path.
% -------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE VARIABLES

file_input = {'YogaTables/yoga1970to1990.tsv', ...
              'YogaTables/yoga1990to2005.tsv', ...
              'YogaTables/yoga2005to2015.tsv', ...
              'YogaTables/yoga2015to2023.tsv'};

% Word used to construct the name of the output directory and to define the
% root of the word dendrogram (generated word tree)
flag_word = upper('YOGA');

% Initialize variables for displaying the progress of processes throughout
% the script.
process_num = 0;
process_num_max = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE PUBMED TABLE AND CREATE PUBMED STRUCT

process_num = process_num + 1;
info = 'PARSE PUBMED TABLE AND CREATE PUBMED STRUCT';
fprintf('[%02i/%02i] %s\n', process_num, process_num_max, info);

PUBMED_STRUCT = parse_pubmed_table(file_input);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD TO PUBMED_STRUCT THE FILTERED WORD LIST

process_num = process_num + 1;
info = 'CREATE FILTERED WORD LIST';
fprintf('[%02i/%02i] %s\n', process_num, process_num_max, info);

PUBMED_STRUCT = create_word_list_filter(PUBMED_STRUCT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE WORD AND TEXT EMBEDDINGS

process_num = process_num + 1;
info = 'CREATE EMBEDDING';
fprintf('[%02i/%02i] %s\n', process_num, process_num_max, info);

PUBMED_STRUCT = create_embbeding(PUBMED_STRUCT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMPORAL CORRELATION ANALYSIS OF WORD OCCURRENCES

process_num = process_num + 1;
info = 'TEMPORAL CORRELATION ANALYSIS OF WORD OCCURRENCES';
fprintf('[%02i/%02i] %s\n', process_num, process_num_max, info);

PUBMED_STRUCT = temporal_correlation(PUBMED_STRUCT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE HTML-TM

process_num = process_num + 1;
info = 'CREATE HTML-TM';
fprintf('[%02i/%02i] %s\n', process_num, process_num_max, info);

output_dir = create_html_tm(PUBMED_STRUCT, flag_word);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE GLOBAL TREE

process_num = process_num + 1;
info = 'CREATE GLOBAL TREE';
fprintf('[%02i/%02i] %s\n', process_num, process_num_max, info);

tree_file = fullfile(output_dir, 'global_tree.nwk');
create_global_tree(PUBMED_STRUCT, flag_word, tree_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
