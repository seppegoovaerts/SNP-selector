%% Setup
clear all; close all; clc; restoredefaultpath;
wd = 'SNP-annotation/';
addpath(genpath([wd 'Functions/']));
cd(wd);
%% 1. Make GENE annotation file
% This file is used used to annotate SNPs to genes.
% 1.1 Load gene list
GeneTable = readtable([wd 'Databases/EnsemblBiomart/NCBI_protCodingGenes_GRCh37_1-23.txt'],'Delimiter','\t');
% 1.2 Match gene names to stalbe IDs
AliasTable = readtable([wd 'Databases/stringdb/9606.protein.aliases.v11.5.txt']);
[ind12,ind21] = vlookupFast(GeneTable.GeneName,AliasTable.alias); 
GeneLegend = table(GeneTable.GeneName(ind21),AliasTable.x_string_protein_id(ind12),'VariableNames',{'Name','stableID'});
GeneLegend.ID = (1:height(GeneLegend))';
% 1.3 Make annotation file
GENES = struct;
GENES.n = length(ind12);
GENES.CHR = GeneTable.Chromosome_scaffoldName(ind21);
GENES.RANGES = [GeneTable.GeneStart_bp_(ind21),GeneTable.GeneEnd_bp_(ind21)];
GENES.NAME = GeneTable.GeneName(ind21);
GENES.buildInfo = 'GRCh37';
% 1.4 Save
save([wd 'Annotation/GENE_ANNOT'],'GENES');
%% 2. Make GO-gene lookup table
% This file will be the main annotation file that can be used with the
% toolbox.
% 2.1 Load Annotations
AnnotationTable = readtable([wd 'Databases//stringdb/9606.protein.enrichment.terms.v11.5.txt']);
% 2.2 Filter for genes in genelist
mask = ismember(AnnotationTable.x_string_protein_id,GeneLegend.stableID);
AnnotationTable = AnnotationTable(mask,:);
[~,idx] = ismember(AnnotationTable.x_string_protein_id,GeneLegend.stableID);
AnnotationTable.geneID = GeneLegend.ID(idx);
% 2.3 Make inventory of gene names and annotations
[unique_annot,idx] = unique(AnnotationTable.term);
unique_genes = unique(GeneLegend.ID);
AnnotationLegend = AnnotationTable(idx,2:4);
AnnotationLegend.descriptionUpper = upper(AnnotationLegend.description);
AnnotationLegend.ID = (1:height(AnnotationLegend))';
[~,idx] = ismember(AnnotationTable.term,AnnotationLegend.term);
AnnotationTable.annotID = AnnotationLegend.ID(idx);
% 2.4 Build matrix
relations = false(length(unique_genes),length(unique_annot));
for i = 1:height(AnnotationTable)
    relations(AnnotationTable.geneID(i),AnnotationTable.annotID(i)) = true;
end
% 2.5 Build dictionaries
gene2annot = dictionary;
for i = 1:size(relations,1)
    gene2annot(i) = {find(relations(i,:))};
end
annot2gene = dictionary;
for i = 1:size(relations,2)
    annot2gene(i) = {find(relations(:,i))};
end
% 2.6 Build struct
ANNOTATION = struct;
ANNOTATION.GENEANNOTATION = GENES;
ANNOTATION.geneInfo = GeneLegend;
ANNOTATION.annotationInfo = AnnotationLegend;
ANNOTATION.annot2gene = annot2gene;
ANNOTATION.gene2annot = gene2annot;
ANNOTATION.nGene = height(GeneLegend);
ANNOTATION.nAnnotation = height(AnnotationLegend);
% 2.7 save
save([wd 'Annotation/ANNOTATION'],'ANNOTATION');
%% END

