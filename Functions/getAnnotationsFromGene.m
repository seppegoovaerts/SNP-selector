function annotations = getAnnotationsFromGene(ANNOTATION,gene)
% ====================================================
% This function finds the annotated terms (GO biological process, tissue
% expression, human phenotype, ...) for each gene in a specified list.
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
% gene - string or cell array of strings containing the (list of) gene(s) to find associated SNPs for
% ====================================================
% OUTPUT:
% ====================================================
% annotations - table with annotations per input gene with table names:
%                    - If more than one input gene is specified, the output is a cell array (1 x nGene), where each element contains a table corresponding to one input gene.
%       - category - cell array (1 x nAnnotation) containing the database / category of each annotation, e.g., 'Biological Process (Gene Ontology)'
%       - term - cell array (1 x nAnnotation) containing the code of each annotation, e.g., 'GO:0019725'
%       - description - cell array (1 x nAnnotation) containing the description of each annotation, e.g., 'Cellular homeostasis'
% ====================================================

%% Funtion body
    % Parse genes
    if ischar(gene), gene = {gene}; end
    nGene = length(gene);
    gene = upper(gene);
    % Get ID
    [~,geneID] = ismember(gene,ANNOTATION.geneInfo.Name);
    % Check invalid gene name
    if nGene == 1 && geneID == 0, warning('Gene not found.'); annotations = []; return; end
    index_not_found = find(geneID == 0);
    for i = 1:length(index_not_found)
        warning(['Could not find: ' gene{index_not_found(i)}]);
    end
    % Lookup
    annotations = cell(1,nGene);
    for i =1:nGene
        if geneID(i) == 0; continue; end
        tmp = ANNOTATION.gene2annot(geneID(i));
        annotationID = tmp{1};
        annotations{i} = ANNOTATION.annotationInfo(annotationID,1:3);
    end
    if nGene == 1, annotations = annotations{1}; end
end
