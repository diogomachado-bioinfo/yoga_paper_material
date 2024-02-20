# WHAT ARE WE LEARNING WITH YOGA: A TEXT MINING APPROACH TO LITERATURE - EXPERIMENT CODE

This repository contains the MATLAB code used to execute the experiment described in the article "What are we learning with Yoga: a text mining approach to literature". The project aims to explore literature in a table with PubMed data through text mining.

## Requirements

- MATLAB R2012a.
- PubMed data table with the following columns (separated by tabs): `pmid` (PubMed identifier), `ti` (title), `ab` (abstract), `fau` (full authors names), and `dp` (date of publication).
- The file available at <https://bit.ly/4bJRazk>. It is a large projection matrix used in the approach. The main script expects this in the current folder or MATLAB path.

## Usage

The `experiment_main_script.m` script manages the execution of the key steps in the experiment.

## Main functions

- `parse_pubmed_table.m`: Primary function to analyze PubMed tables, expecting columns `pmid`, `title`, `abstract`, and `authors` (separated by tabs).
- `create_word_list_filter.m`: Main function to filter words based on specific criteria.
- `create_embbeding.m`: Key function to create word and text embeddings.
- `temporal_correlation.m`: Primary function to analyze the temporal correlation of word occurrences.
- `create_html_tm.m`: Main function to generate HTML files based on PubMed data.
- `create_global_tree.m`: Main function to build a global dendrogram of words (word tree).

## Other functions

There are additional functions used in this project. Refer to the individual function files for details.

## Authors

- Roberto T. Raittz (Code Originality)
- Diogo de J. S. Machado (Code Organization)
