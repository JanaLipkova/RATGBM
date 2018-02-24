function SyntheticCaseParametricIC

close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/')

% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
bSynthetic = 1;

GenId = 14;
myfilename = sprintf('../../../TMCMC/SmallSyntheticParametricIC/SyntheticAll/SynthethicAll_PIC_12K/curgen_db_%03d.txt',GenId);
mydata = importdata(myfilename);
[Nx,Ny] = size(mydata);

names =       ['D   ' ;'rho '  ;'Tend' ;'r   '; 'alph'; 'beta';'PETn';'b   ';'T1uc';'T2uc';'Tn  '];
groundTruth = [1.3e-03, 2.5e-02, 302, 0.0062, 312.55,   30.43,  0.019,0.87, 0.7,    0.25,  5.0e-02];


for i = 1:Ny-1
    mydata(:,i) = exp( mydata(:,i) );
end;

B1	=   exp([	-8.9480   -3.2702 ]);
B2	=   exp([	-5.9145   -1.6607 ]);
B3	=   exp([	 3.4012    7.3132 ]);
B4  =   exp([   -6.9077   -3.1881 ]);
B5  =   exp([   -6.9077    5.8862 ]);
B6  =   exp([   -6.9077    5.1930 ]);
B7	=   exp([	-4.6052   -0.9163 ]);
B8	=   exp([	-0.9163    0.0488 ]);
B9	=   exp([	-0.5108   -0.2231 ]);
B10	=   exp([	-4.6052   -0.9163 ]);
B11	=   exp([	-2.9957   -2.3026 ]);

bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9;B10;B11];

fname = sprintf('S_pIC_gen_%d_%i.txt',GenId,Nx);
fid = fopen(fname, 'wt'); % Open for writing
[Nx,Ny] = size(mydata);
for i=1:Nx
    for j=1:Ny
        fprintf(fid, '%d ', mydata(i,j));
    end
    fprintf(fid, '\n ');
end;
fclose(fid);



%%
% bestC = find( max(mydata(:,end)) == mydata(:,end));
% best = mydata(bestC,:)
% meanData = mean(mydata)
% varData = var(mydata)
% stdData = sqrt(varData)


% 4) Plot results - Marginals, KDE and scattered samples
%----------------------------------------------------------------
KDEbins = 15*ones(1,Ny-1);
Nbins   = 200;
dump = 1;
% param= [1,2,3,4,5,6];%1:Ny-2 ;
param= [1,7,8,9,10,11];%1:Ny-2 ;


plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)
