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


% function TMCMC_evaluate_paramterSpace
close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/')


% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
bSynthetic = 0;
Patient = 1;


GenId = 9;

myfilename = sprintf('../../../TMCMC/PatientCases/Patient01/Patient01PET/curgen_db_%03d.txt',GenId);
mydata = importdata(myfilename);

%rescale
mydata(:,1) = exp( mydata(:,1).* 100 );
mydata(:,2) =      mydata(:,2).* 10;
mydata(:,3) = exp( mydata(:,3).* 100 );
mydata(:,4) =      mydata(:,4)       ;
mydata(:,5) = exp( mydata(:,5)* 100 );
mydata(:,6) =      mydata(:,6)*100;




fid = fopen('P01_curgen_db_09.txt', 'wt'); % Open for writing
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



if (bSynthetic)
    
    B1       = [   -9.0e-02   -3.2e-02] * 100;
    B2       = [    2.7e-04    1.9e-02] * 10;
    B3       = [    4.1e-02    7.0e-02] * 100;   % Synthetic
    B4       = [    1.0e-04,   3.7e-02] * 10;    % Synthetic
    B5       = [   -4.6e-02   -1.2e-02] * 100;
    B6       = [    4.0e-03    1.0e-02] * 100;
    
    
else 
    if (Patient == 1)
        B1       = [   -9.0e-02   -3.2e-02] * 100;
        B2       = [    2.7e-04    1.9e-02] * 10;
        B3       = [    3.4e-02    7.2e-02] * 100;   % P01
        B4       = [    3.0e-03    2.0e-02]      ;
        B5       = [   -4.6e-02   -1.2e-02] * 100;
        B6       = [    4.0e-03    1.0e-02] * 100;
    else
        B1       = [   -9.0e-02   -3.2e-02] * 100;
        B2       = [    2.7e-04    1.9e-02] * 10;
        B3       = [    4.2e-02    7.7e-02 ] * 100; %P11
        B4       = [    1.0e-03    8.0e-02]      ;
        B5       = [   -4.6e-02   -1.2e-02] * 100;
        B6       = [    4.0e-03    1.0e-02] * 100;
    end;
end;
B1       = exp(B1);
B3       = exp(B3);
B5       = exp(B5);

bounds = [B1;B2;B3;B4;B5;B6];



% 4) Plot results - Marginals, KDE and scattered samples
%----------------------------------------------------------------
 KDEbins = 15*ones(1,6);
%  KDEbins = [15, 15, 15, 10, 15, 20];  % P11
% % KDEbins = [10, 10, 20, 10, 50, 50];  % P11
% KDEbins = [10, 10, 20, 10, 20, 10];  % P01

Nbins   = 200;
dump = 1;
param=[1,2, 3,4,5,6];

plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)
