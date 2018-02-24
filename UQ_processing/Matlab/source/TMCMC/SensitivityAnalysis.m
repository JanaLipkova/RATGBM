%=====================================================
%
%   Sensitivity analysis of model parameters
%
%------------------------------------------------------
% SENSITIVE PARAMETERS
% - parameters are inferred quickly (i.e in earlier generations)
% - and have narrow credible intervals
% LESS SENSITIVE PARAMETERS
% - parameters get inferred in later populations
% - and are not very localised by the posterior distribution
% NOT INFERRABLE PARAMETER GIVEN THE AVAILABLE DATA:
% - if the distribution doesn?t change much between generations and resemble uniform distr. from 0th generation
% =====================================================

% function SensitivityAnalysis

close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/')



GenId = [4,8,12,16];
pathToData = '../../../TMCMC/PatientCases/Patient01/SensitivityTestAll/Patient01All_8K/';
names =       ['D   ' ;'rho '  ;'Tend' ;'ICn '; 'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
scaling = [1, 10, 100000, 1, 100, 100, 100, 100, 10];


B1  = exp ([ 	-8.9480   -3.2702] ) .*scaling(1);
B2  = exp ([	-8.2171   -3.9633] ) .*scaling(2);
B3	= exp ([	-8.1117   -4.2687] ) .*scaling(3);
B4	= exp ([	-5.8091   -3.9120] ) .*scaling(4);
B5	= exp ([	-9.2103   -5.8091] ) .*scaling(5);
B6	= exp ([	-5.5215   -4.6052] ) .*scaling(6);
B7	= exp ([	-5.1160   -4.8283] ) .*scaling(7);
B8	= exp ([	-9.2103   -5.5215] ) .*scaling(8);
B9	= exp ([	-5.2983   -4.6052] ) .*scaling(9);

bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9];
param= 1:length(scaling);

plotSensitivityAnalysis(pathToData, GenId, param, names, bounds, scaling )
