%======================================================
%
%  Function to plot currgen_db files from TMCMC
%
%======================================================
%
% INPUT:
%   1) range of parameter values (as defined in tmcmc.par)
%   2) Path to the curgen_db file
%======================================================

% function TMCMC_evalutate_current_gen

close all; clear all; clc;

addpath('../../lib/jbfill')
addpath('../../lib/kde2d')
addpath('../../lib/')


% 1) range of parameter values
% 2) Names of parameres & Ground Trurh

smallerRange = 0;
bSynthetic = 1;  % 1 for synthetic data i.e also ground truth, 0 for patient data

if(smallerRange == 0)
    
%     B1 =[ 1.3e-04   3.8e-02 ];
    B1 =[-8.9e-02   -3.2e-02];
    B2 =[ 2.7e-04	1.9e-02 ] * 10;
    B3 =[ 3.0e-04   1.05e-02] * 100000;
    B4 =[ 1.0e-03   3.7e-02 ] * 10;
    B5 =[ -4.6e-02 -1.2e-02 ] * 100;
    B6 =[ 5.0e-03   1.0e-02 ] * 100;
    
    names =       ['D   ' ;'rho '  ;'Tend'   ;'ICn ' ;'PETn';'b   '    ];
%     groundTruth = [0.0013, 0.025, 270, 0.05, log(0.05) , 0.827702  ];  % GT for Bigger Added Noise
    groundTruth = [log(0.0013), 0.025, 270, 0.05, log(0.0251) , 0.6793  ];  % GT for Bigger Added Noise

else
    
%     B1 =[ 1.3e-04   3.8e-02 ];
    B1      = [1.3e-04   1.3e-02];
    B2      = [2.7e-04	  1.9e-02]  *10;
    
    B3      = [6.0e-04    1.05e-02] * 100000;
    
%     B3      = [9.0e-04   5.0e-03]   * 100000 ;
    B4      = [1.0e-03   3.7e-02]  * 10;
    B5      = [-4.6e-02   -1.2e-02] * 100;
    B6      = [ 5.0e-03   1.0e-02] * 100;
    
    names =       ['D   ' ;'rho '  ;'Tend'   ;'ICn ' ;'PETn';'b   '    ];
      groundTruth = [0.0013, 0.025, 270, 0.01, log(0.0251) , 0.6793  ];  % GT for Smaller Added Noi
%      groundTruth = [0.0013, 0.025, 270, 0.01, 0.02, 0.6793  ];  % GT for Smaller Added Noise

end

bounds = [B1;B2;B3;B4;B5;B6];




%3) Plot statistic for pairs of parameters on given Path
parameters = [5,6];
dump = 1;

KDEbins = 20*ones(1,6);
Nbins   = 200;

GenId =10;
myfilename = sprintf('../../../TMCMC/SynthethicPETdataSmallerAddedNoiseLogD/curgen_db_%03d.txt',GenId);

mydata = importdata(myfilename);


mydata(:,2) = mydata(:,2).*10;
mydata(:,3) = mydata(:,3).*100000;
mydata(:,4) = mydata(:,4).*10;
mydata(:,5) = mydata(:,5)*100;
mydata(:,6) = mydata(:,6)*100;

plotTMCMC_PairStatistic(mydata,parameters,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)
%
%
% %% coordinates of the best estimate:
bestC = find( max(mydata(:,7)) == mydata(:,7));
best = mydata(bestC,:)
meanData = mean(mydata)
varData = var(mydata)
stdData = sqrt(varData)


% close all;
% figure,
% data = [mydata(:,parameters(1)), mydata(:,parameters(2))];
% MIN_XY = B1;
% MAX_XY = B2;
%
%
%  [bandwidth,density,X,Y]=kde2d(data);
%  % plot the data and the density estimate
%  surf(X,Y,density,'LineStyle','none')
%  colormap hot, hold on, alpha(.8)
%  set(gca, 'color', 'blue');
%  plot(data(:,1),data(:,2),'w.','MarkerSize',5)

% [bandwidth,density,X,Y]=kde2d(data,2^8,MIN_XY,MAX_XY);
% pcolor(X,Y,density);
% shading flat
% figure,
%
% data
% % apply routine
%  [bandwidth,density,X,Y]=kde2d(data);
%  % plot the data and the density estimate
%  surf(X,Y,density,'LineStyle','none'), view([0,70])
%  colormap hot, hold on, alpha(.8)
%  set(gca, 'color', 'blue');
%  plot(data(:,1),data(:,2),'w.','MarkerSize',5)
