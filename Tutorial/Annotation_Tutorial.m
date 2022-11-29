%% Setup
clear all; close all; clc; restoredefaultpath;
wd = 'SNP-annotation/';
cd(wd);
%% Load toolbox functions and annotation file
addpath(genpath([wd 'Functions/']));
load([wd 'Annotation/ANNOTATION.mat']);
%% 1. Get genes from annotation
% You can look-up by both description (case insensitive) and code.
% Example 1: "GO Biological processes"
genes = getGenesFromAnnotation(ANNOTATION,{'intramembranous ossification','GO:0010763'});

% Example 2: "Tissue expression"
genes = getGenesFromAnnotation(ANNOTATION,'bone');

%% 2. Get annotation from genes
% Example 1:
annotations = getAnnotationsFromGene(ANNOTATION,{'SOX9','RUNX2'});

%% 3. Get SNPS from genes
% Load Heel Bone Mineral Density GWAS summary
sumstats = readtable([wd 'Tutorial/GWAS_summary_Heel_Bone_Mineral_Density_UKBB.txt']);

% Example 1:
% Make gene list
genelist = {'BMP2','EYA4','DLX5','DLX6','RUNX2','WNT16','EN1','LRP5','BMP2','BMPR2','KMT2D','TBX1','RSPO3','EMP1','ESR1','CPED1'};
% Find SNPs
margin = 500e3;
index = getSnpsFromGene(ANNOTATION.GENEANNOTATION,sumstats.chr,sumstats.pos,genelist,margin);
% Visualize on Manhattan plot
genomewideplot(sumstats.chr,sumstats.pos,sumstats.logP,5e-8,index)

%% 4. Get SNPS from annotations
% Example 1 "Tissue expression: bone":
% Make annotation list
annotationlist = {'bone'};
% Set gene margin
margin = 500e3;
% Gind SNPs
index = getSnpsFromAnnotation(ANNOTATION,sumstats.chr,sumstats.pos,annotationlist,margin);
% Visualize on Manhattan plot
genomewideplot(sumstats.chr,sumstats.pos,sumstats.logP,5e-8,index)

% Example 2: "A single GO Biological process: skeletal system development"
% Make annotation list
annotationlist = {'skeletal system development'};
% Set gene margin
margin = 500e3;
% Gind SNPs
index = getSnpsFromAnnotation(ANNOTATION,sumstats.chr,sumstats.pos,annotationlist,margin);
% Visualize on Manhattan plot
genomewideplot(sumstats.chr,sumstats.pos,sumstats.logP,5e-8,index)

% Example 3: "GO Biological processes related to ossification"
% Make annotation list
annotationlist = {'ossification','osteoblast differentiation','intramembranous ossification','regulation of ossification','ossification involved in bone remodeling'};
% Set gene margin
margin = 500e3;
% Gind SNPs
index = getSnpsFromAnnotation(ANNOTATION,sumstats.chr,sumstats.pos,annotationlist,margin);
% Visualize on Manhattan plot
genomewideplot(sumstats.chr,sumstats.pos,sumstats.logP,5e-8,index)

% Example 4: "Human phenotype: Abnormal long bone morphology"
% Make annotation list
annotationlist = {'Abnormality of long bone morphology'};
% Set gene margin
margin = 500e3;
% Gind SNPs
index = getSnpsFromAnnotation(ANNOTATION,sumstats.chr,sumstats.pos,annotationlist,margin);
% Visualize on Manhattan plot
genomewideplot(sumstats.chr,sumstats.pos,sumstats.logP,5e-8,index)

% Example 5: "Human phenotype: Osteolysis"
% Make annotation list
annotationlist = {'osteolysis'};
% Set gene margin
margin = 500e3;
% Gind SNPs
index = getSnpsFromAnnotation(ANNOTATION,sumstats.chr,sumstats.pos,annotationlist,margin);
% Visualize on Manhattan plot
genomewideplot(sumstats.chr,sumstats.pos,sumstats.logP,5e-8,index)

%% The End
