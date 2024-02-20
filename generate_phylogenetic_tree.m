function [phylogenetic_tree, distance_matrix] = generate_phylogenetic_tree(data_matrix, names, varargin)
% Generates a phylogenetic tree from the input data matrix using hierarchical clustering.
%
% Input:
%   - data_matrix: Matrix containing the data for phylogenetic analysis.
%   - names: Names of the elements in the matrix.
%   - varargin:
%       - 'Distance': Distance metric for clustering (default: 'Euclidean').
%       - 'Method': Clustering method (default: 'ward').
%       - 'Pdist': Precomputed distance matrix (optional).
%       - 'Show': Flag to show the tree visualization (default: 1).
%       - 'Reroot': Node to be used as the root in the tree (optional).
%
% Output:
%   - phylogenetic_tree: Phylogenetic tree structure.
%   - distance_matrix: Computed distance matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = varargin;
distance_metric = find_in_varargin(V, 'Distance', 'Euclidean');
clustering_method = find_in_varargin(V, 'Method', 'ward');
precomputed_distance = find_in_varargin(V, 'Pdist', []);
show_tree = find_in_varargin(V, 'Show', 1);
reroot_node = find_in_varargin(V, 'Reroot', []);

if ~isempty(precomputed_distance)
    distance_matrix = precomputed_distance;
else
    if isnumeric(distance_metric)
        custom_distance = eval(['@(x, y) Lkdist(x, y,' num2str(distance_metric) ')']);
        distance_matrix = pdist(data_matrix, custom_distance);
    else
        distance_matrix = pdist(data_matrix, distance_metric);
    end
end

if ~strcmpi(clustering_method, 'NJ')
    linkage_result = linkage(distance_matrix, clustering_method);
    phylogenetic_tree = phytree(linkage_result, names);
else
    phylogenetic_tree = seqneighjoin(distance_matrix, 'equivar', names);
end

if ~isempty(reroot_node)
    if iscell(reroot_node)
        node_to_reroot = getbyname(phylogenetic_tree, reroot_node{1});
    else
        node_to_reroot = getbyname(phylogenetic_tree, reroot_node);
    end
    
    if sum(node_to_reroot) > 1
        take_one = find(node_to_reroot, 1);
        node_to_reroot = node_to_reroot ~= node_to_reroot;
        node_to_reroot(take_one) = 1;
    end
    
    if sum(node_to_reroot) == 1
        phylogenetic_tree = reroot(phylogenetic_tree, node_to_reroot);
    end
end

if show_tree == 1
    phytreetool(phylogenetic_tree);
end

end

function distance_matrix = Lkdist(XI, XJ, varargin)
% Generalized Lk distance function with parameter k for Euclidean distance.
%
% Input:
%   - XI, XJ: Input vectors or matrices for distance calculation.
%   - varargin:
%       - k: Exponent for the distance calculation (default: 2).
%
% Output:
%   - distance_matrix: Computed Lk distance matrix.

V = varargin;
if ~isempty(V)
    k = V{1};
else
    k = 2;
end

squared_diff = abs(repmat(XI, length(XJ(:,1)), 1) - XJ).^k;
D2_squared = sum(squared_diff, 2);
distance_matrix = (D2_squared).^(1/k);

end

