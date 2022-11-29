function [index,chr_filtered,pos_filtered,gene_names] = getSnpsFromGene(GENEANNOTATION,chr,pos,gene,margin)
% ====================================================
% This function finds the index of SNPs within a specified margin of the
% specified list of genes.
% ====================================================
% INPUT:
% ====================================================
% GENEANNOTATION - struct object with fields:
%       - NAME: cell array (1 x nGene) containing one name (string) per gene
%       - CHR: numeric array (1 x nGene) containing the chromosome number of each gene
%       - RANGES: numeric array (2 x nGene) containing start and end position of each gene
%       -n: integer value indicating the total number of genes in the GENEANNOTATION object
% chr - numeric array (1 x nSNP) containing the chromosome number of each input SNP
% pos - numberic array (1 x nSNP) containing the chromosomal position of each input SNP
% gene - string or cell array of strings containing the (list of) gene(s) to find associated SNPs for
% margin - numeric value indicating the maximal distance (in base pairs) to a gene a SNP can be (default: 500 kb)
% ====================================================
% OUTPUT:
% ====================================================
% index - numeric array (1 x n_classified_snps) containing the indices of SNPs classified to the specified genes
% chr_filtered - numeric array (1 x n_classified_snps) containing the chromosome numbers of SNPs classified to the specified genes
% pos_filtered - numeric array (1 x n_classified_snps) containing the chromosomal positions of SNPs classified to the specified genes
% gene_names - cell array (1 x n_classified_snps) containing the gene each
% associated to each classified SNP (warning: this can be overwritten, causing the final annotation to not necessarily reflect the closest gene for SNPs that fall within the specified margin of two genes)
% ====================================================

%% Funtion body
if nargin < 4, margin = 500e3; end
    % Parse input
    if numel(chr) ~= numel(pos), error('Chromosome number and position vectors must have equal length.'); end
    if ischar(gene), gene = {gene}; end
    nSNP = length(chr);
    nGene = length(gene);
    gene = upper(gene);
    % Check invalid gene name
    [~,geneID] = ismember(gene,GENEANNOTATION.NAME);
    index_not_found = find(geneID == 0);
    for i = 1:length(index_not_found)
        warning(['Could not find: ' gene{index_not_found(i)}]);
    end
    GA = filterGeneAnnotation(GENEANNOTATION,gene);
    % Find SNPs
    found_mask = false(nSNP,1);
    gene_idx = zeros(nSNP,1);
    for i = 1:GA.n
        in_gene_mask = (chr == GA.CHR(i)) & (pos > GA.RANGES(i,1) - margin) & (pos < GA.RANGES(i,2) + margin);
        found_mask = found_mask | in_gene_mask;
        gene_idx(in_gene_mask) = i;
    end
    % Output
    index = find(found_mask);
    if nargout > 1
        chr_filtered = chr(index);
        pos_filtered = pos(index);
        gene_names = GA.NAME(gene_idx(found_mask));
    end
end

%% Helper
function SUBGENE = filterGeneAnnotation(ANNOT,GENELIST)
    m = ismember(ANNOT.NAME,GENELIST);
    if sum(m) == 0, error('Gene(s) not found in annotation file.'); end
    SUBGENE.CHR = ANNOT.CHR(m);
    SUBGENE.RANGES = ANNOT.RANGES(m,:);
    SUBGENE.NAME = ANNOT.NAME(m);
    SUBGENE.n = sum(m);
end