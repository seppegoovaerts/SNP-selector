function [index,chr_filtered,pos_filtered,gene_names] = getSnpsFromAnnotation(ANNOTATION,chr,pos,annotation,margin)
% ====================================================
% This function finds the index of SNPs within a specified margin of the
% the genes annotated with one of the specified terms (GO biological
% process, tissue expression, human phenotype, ...).
% ====================================================
% INPUT:
% ====================================================
% ANNOTATION - struct object with fields:
%       - geneInfo - "Database of genes" - table with names:
%               - Name - cell array (1 x nGene) containing one name (string) per gene
%               - stableID - cell array (1 x nGene) containing the StringDB stable protein ID (string) corresponding to each gene
%               - ID - numeric array (1 x nGene) containing an integer ID for each gene for easy look-ups
%       - annotationInfo - ''Database of annotations' - table with names:
%               - category - cell array (1 x nAnnotation) containing the database / category of each annotation, e.g., 'Biological Process (Gene Ontology)'
%               - term - cell array (1 x nAnnotation) containing the code of each annotation, e.g., 'GO:0019725'
%               - description - cell array (1 x nAnnotation) containing the description of each annotation, e.g., 'Cellular homeostasis'
%               - descriptionUpper - cell array (1 x nAnnotation) containing the description of each annotation in uppercase for easier look-ups
%               - ID - numeric array (1 x nAnnotation) containing an integer ID for each annotation for easy look-ups
%       - annot2gene - "Gene look-up table" - dictionary:
%               - keys: annotation ID (ANNOTATION.annotationInfo.ID)
%               - values: gene ID (ANNOTATION.geneInfo.ID)
%       - gene2annot - "Annotation look-up table" - dictionary:
%               - keys: gene ID (ANNOTATION.geneInfo.ID)
%               - values: annotation ID (ANNOTATION.annotationInfo.ID)
%       - GENEANNOTATION - struct object with fields:
%               - NAME: cell array (1 x nGene) containing one name (string) per gene
%               - CHR: numeric array (1 x nGene) containing the chromosome number of each gene
%               - RANGES: numeric array (2 x nGene) containing start and end position of each gene
%               -n: integer value indicating the total number of genes in the GENEANNOTATION object
% chr - numeric array (1 x nSNP) containing the chromosome number of each input SNP
% pos - numberic array (1 x nSNP) containing the chromosomal position of each input SNP
% annotation - string or cell array of strings containing the (list of) annotations(s) to find associated SNPs for. These can either be codes or descriptions.
% margin - numeric value indicating the maximal distance (in base pairs) to a gene a SNP can be (default: 500 kb)
% ====================================================
% OUTPUT:
% ====================================================
% index - numeric array (1 x n_classified_snps) containing the indices of SNPs classified to the specified genes
% chr_filtered - numeric array (1 x n_classified_snps) containing the chromosome numbers of SNPs classified to the specified genes
% pos_filtered - numeric array (1 x n_classified_snps) containing the chromosomal positions of SNPs classified to the specified genes
% gene_names - cell array (1 x n_classified_snps) containing the gene each associated to each classified SNP (warning: this can be overwritten, causing the final annotation to not necessarily reflect the closest gene for SNPs that fall within the specified margin of two genes)
% ====================================================

%% Funtion body
if nargin < 4, margin = 500e3; end
    % Parse input
    if numel(chr) ~= numel(pos), error('Chromosome number and position vectors must have equal length.'); end
    if ischar(annotation), annotation = {annotation}; end
    % Get genes
    genes = getGenesFromAnnotation(ANNOTATION,annotation);
    if isempty(genes), error('No genes found for given annotations.'); end
    % Combine genes
    if iscell(genes{1})
        all_genes = cat(1,genes{:});
    else
        all_genes = genes;
    end
    % Get SNPs
    [index,chr_filtered,pos_filtered,gene_names] = getSnpsFromGene(ANNOTATION.GENEANNOTATION,chr,pos,all_genes,margin);
end

