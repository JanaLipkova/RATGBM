%==========================================
%
%  read in samples from all generations,
%  rescale them to correct parameter range
%  and put all in one file
%
%==========================================


% function combineAllsamples

close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/')

j = 1;
for GenId =7:8
    
    myfilename = sprintf('../../../TMCMC/PatientCases/Patient01/Patient01PET/curgen_db_%03d.txt',GenId);
    mydata = importdata(myfilename);
    
    %rescale
    mydata(:,1) = exp( mydata(:,1).* 100 );
    mydata(:,2) =      mydata(:,2).* 10;
    mydata(:,3) = exp( mydata(:,3).* 100 );
    mydata(:,4) =      mydata(:,4)       ;
    mydata(:,5) = exp( mydata(:,5)* 100 );
    mydata(:,6) =      mydata(:,6)*100;
    
    [Nx,Ny] = size(mydata);
    alldata( Nx*(j -1 )+ 1 : Nx*j, :) = mydata;
    j=j+1;
end;



fid = fopen('P01_all.txt', 'wt'); % Open for writing
[Nx,Ny] = size(alldata);
for i=1:Nx
    for j=1:Ny
        fprintf(fid, '%d ', alldata(i,j));
    end
    fprintf(fid, '\n ');
end;
fclose(fid);

%%
names =       ['D   ' ;'rho '  ;'Tend'   ;'ICn ' ;'PETn';'b   '];

B1       = [   -9.0e-02   -3.2e-02] * 100;
B2       = [    2.7e-04    1.9e-02] * 10;
B3       = [    3.4e-02    7.2e-02] * 100;   % P01
B4       = [    3.0e-03    2.0e-02]      ;
B5       = [   -4.6e-02   -1.2e-02] * 100;
B6       = [    4.0e-03    1.0e-02] * 100;

B1       = exp(B1);
B3       = exp(B3);
B5       = exp(B5);

bounds = [B1;B2;B3;B4;B5;B6];

% 4) Plot results - Marginals, KDE and scattered samples
%----------------------------------------------------------------
 KDEbins = 15*ones(1,6);

Nbins   = 200;
dump = 1;
param=[1,2, 3,4,5,6];
bSynthetic=0;
GenId=8;
groundTruth=0;

plotTMCMC_AllStatistic(alldata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)


