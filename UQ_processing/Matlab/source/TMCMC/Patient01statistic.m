%function Patient01statistic

% function TMCMC_evaluate_paramterSpace
close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/')

% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
dataType = 'All';

if(dataType == 'PET')
    
    GenId = 9;
    myfilename = sprintf('../../../TMCMC/PatientCases/Patient01/SensitivityTestPET/Patient01PET_10K/curgen_db_%03d.txt',GenId);
        
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend'   ;'ICn ' ;'PETn';'b   '];
    
    %rescale
    scaling = [1, 10, 100000, 1, 100, 100];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) ) .* scaling(i);
    end;

    
    % range
    B1 	= exp( [	-8.9480   -3.2702 ] ) .*scaling(1);
    B2	= exp( [	-8.2171   -3.9633 ] ) .*scaling(2);
    B3	= exp( [	-8.1117   -4.2687 ] ) .*scaling(3);
    B4	= exp( [	-5.8091   -3.9120 ] ) .*scaling(4);
    B5	= exp( [	-9.2103   -5.8091 ] ) .*scaling(5);
    B6	= exp( [	-5.5215   -4.6052 ] ) .*scaling(6);
    
    bounds = [B1;B2;B3;B4;B5;B6];
    
    fname = sprintf('P01_PET_gen_%d_%i.txt',GenId,Nx);
    fid = fopen(fname, 'wt'); % Open for writing
    [Nx,Ny] = size(mydata);
    for i=1:Nx
        for j=1:Ny
            fprintf(fid, '%d ', mydata(i,j));
        end
        fprintf(fid, '\n ');
    end;
    fclose(fid);
end;


if(dataType == 'All')
    
    GenId = 17;
    myfilename = sprintf('../../../TMCMC/PatientCases/Patient01/SensitivityTestAll/Patient01All_12K/curgen_db_%03d.txt',GenId);
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'ICn '; 'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    
    %rescale
    scaling = [1, 10, 100000, 1, 100, 100, 100, 100, 10];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) ) .* scaling(i);
    end;

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
    
    fname = sprintf('P01_ALL_gen_%d_%i.txt',GenId,Nx);
    fid = fopen(fname, 'wt'); % Open for writing
    [Nx,Ny] = size(mydata);
    for i=1:Nx
        for j=1:Ny
            fprintf(fid, '%d ', mydata(i,j));
        end
        fprintf(fid, '\n ');
    end;
    fclose(fid);
end;


%%
bestC = find( max(mydata(:,end)) == mydata(:,end));
best = mydata(bestC,:)
meanData = mean(mydata)
varData = var(mydata)
stdData = sqrt(varData)


% 4) Plot results - Marginals, KDE and scattered samples
%----------------------------------------------------------------
KDEbins = 15*ones(1,Ny-1);
Nbins   = 200;
dump = 1;
param= 1:Ny-1 ;

bSynthetic = 0;
groundTruth = 0;

plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)
