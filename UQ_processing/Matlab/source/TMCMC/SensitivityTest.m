%======================================================
%
%  For given curgen_db file from TMCMC create join plot including:
%  1) marginals (diagonal)
%  2) scatter plot for samlpes (above diagonal)
%  3) Gaussian KDE of samples (belov diagonal)
%
%======================================================
%
% INPUT:
%   1) Path to the curgen_db file + Index of generation to be plot
%   2) range of parameter values (as define in tmcmc.par)
%   3) If synthetic data - ground truth
%======================================================

% function SensitivityTest

close all; clear all; clc
addpath('../../lib/jbfill')
addpath('../../lib/')


% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
bSynthetic = 1;
GenId = 11;
myfilename = sprintf('../../../TMCMC/SensitivityTest/SynthethicPET_allLog_5K/curgen_db_%03d.txt',GenId);
mydata = importdata(myfilename);

%rescale
mydata(:,1) = exp( mydata(:,1) )    ;
mydata(:,2) = exp( mydata(:,2) ).* 10;
mydata(:,3) = exp( mydata(:,3) ).* 100000;
mydata(:,4) = exp( mydata(:,4) ).* 10;
mydata(:,5) = exp( mydata(:,5) ).* 100;
mydata(:,6) = exp( mydata(:,6) ).* 100;

fid = fopen('S_5K_curgen_db_11.txt', 'wt'); % Open for writing
[Nx,Ny] = size(mydata);
for i=1:Nx
    for j=1:Ny
        fprintf(fid, '%d ', mydata(i,j));
    end
    fprintf(fid, '\n ');
end;
fclose(fid);


%% 3) Range of paramert values + GroundTruth + Names of parameters
%----------------------------------------------------------------
names =       ['D   ' ;'rho '  ;'Tend'   ;'ICn ' ;'PETn';'b   '];
groundTruth = [0.0013, 0.025, 270, 0.0058, 0.0251 , 0.6793  ];

B1 =[ -8.9479     -3.2701 ];
B2 =[ -8.2170     -3.9633 ];
B3 =[ -8.1117     -4.5098 ];
B4 =[ -9.2103     -3.2968 ];
B5 =[ -9.2103	  -5.8091 ];
B6 =[ -5.5214     -4.6051 ];

B1 = exp(B1);
B2 = exp(B2).* 10;
B3 = exp(B3).* 100000;
B4 = exp(B4).* 10;
B5 = exp(B5).* 100;
B6 = exp(B6).* 100;

bounds = [B1; B2; B3; B4; B5; B6];



% 4) Plot results - Marginals, KDE and scattered samples
%----------------------------------------------------------------
KDEbins = 15*ones(1,6);
Nbins   = 200;
dump    = 1;
param   = 1:Ny-1;

plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)

