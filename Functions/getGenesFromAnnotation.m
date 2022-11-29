function genes = getGenesFromAnnotation(ANNOTATION,annotation)
% ====================================================
% This function finds the genes with specified annotations (GO
% biological process, tissue expression, human phenotype, ...) .
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
% annotation - string or cell array of strings containing the (list of) annotations(s) to find associated SNPs for. These can either be codes or descriptions.
% ====================================================
% OUTPUT:
% ====================================================
% genes - cell array (1 x nGene) with gene names
%                    - If more than one input annotation is specified, the output is a cell array (1 x nAnnotation), where each element
%                      contains a cell array with gene names for that annotation.
% ====================================================

%% Funtion body
    % Parse annotations
    if ischar(annotation), annotation = {annotation}; end
    nAnnot = length(annotation);
    annotation = upper(annotation);
    % Check for duplicate annotation descriptions
    [annotation_matched_mask,input_index] = ismember(ANNOTATION.annotationInfo.descriptionUpper,annotation);
    input_index_found = input_index(annotation_matched_mask);
    [gc,gr] = groupcounts(input_index_found);
    duplicate_input_index = gr(gc > 1);
    if any(duplicate_input_index)
        warning('One or more descriptions match multiple annotations. Consider searching by code instead.');
        disp('Overview of duplicate matches:')
        duplicate_annotation_mask = ismember(input_index,duplicate_input_index);
        disp(ANNOTATION.annotationInfo(duplicate_annotation_mask,[2,1,3]))
    end
    % Get ID
    [~,annotationID1] = ismember(annotation,ANNOTATION.annotationInfo.descriptionUpper);
    [~,annotationID2] = ismember(annotation,ANNOTATION.annotationInfo.term);
    annotationID = annotationID1 + annotationID2;
    % Check for unmatched processes
    if nAnnot == 1 && annotationID == 0, warning('Annotation not found.'); genes = []; return; end
    unmatched_mask = annotationID == 0;
    if any(unmatched_mask)
        warning('One or more annotation terms could not me matched to database.');
        disp('Could not find the following terms:')
        fprintf([strjoin(annotation(unmatched_mask),'\n'),'\n'])
    end
    % Lookup
    genes = cell(1,nAnnot);
    for i =1:nAnnot
        if annotationID(i) == 0; continue; end
        tmp = ANNOTATION.annot2gene(annotationID(i));
        geneID = tmp{1};
        genes{i} = ANNOTATION.geneInfo.Name(geneID);
    end
    if nAnnot == 1, genes = genes{1}; end
end
