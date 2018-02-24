% function SyntheticSmallCase


close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/')

% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
dataType = 'All';
bOutputR = 0;

if(dataType == 'PET')
    
    GenId = 14;
    myfilename = sprintf('../../../TMCMC/SmallSynthetic/SynthethicPETonly/SynthethicPETonly_8K/curgen_db_%03d.txt',GenId);
        
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);    
    names =       ['D   ' ;'rho '  ;'Tend'   ;'ICn ' ;'PETn';'b   '];
    groundTruth = [1.3e-03, 2.5e-02, 302, 0.0062,    0.019,  0.87];

    %rescale
    scaling = [1, 1, 1, 1, 1, 1];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) ) .* scaling(i);
    end;
    
    % range
    B1   =   exp([   -8.9480   -3.2702 ]);
    B2   =   exp([   -5.9145   -1.6607 ]);
    B3   =   exp([    3.4012    7.3132 ]);
    B4   =   exp([   -6.9078   -1.2040 ]);
    B5   =   exp([   -4.6052   -0.9163 ]);
    B6   =   exp([   -0.9163    0.0488 ]);
    
    bounds = [B1;B2;B3;B4;B5;B6];
    
    if(bOutputR)
        fname = sprintf('SS_PET_gen_%d_%i.txt',GenId,Nx);
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

end;

if(dataType == 'All')
    
    GenId = 14;
    myfilename = sprintf('../../../TMCMC/SmallSynthetic/SynthethicAll/SynthethicAll_14K/curgen_db_%03d.txt',GenId);
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'ICn '; 'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    groundTruth = [1.3e-03, 2.5e-02, 302, 0.0062,    0.019, 0.87, 0.7,     0.25,  5.0e-02];

    %rescale
    scaling = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) ) .* scaling(i);
    end;
    
    B1   =   exp([   -8.9480   -3.2702 ]);
    B2   =   exp([   -5.9145   -1.6607 ]);
    B3   =   exp([    3.4012    7.3132 ]);
    B4   =   exp([   -6.9078   -1.2040 ]);
    B5   =   exp([   -4.6052   -0.9163 ]);
    B6   =   exp([   -0.9163    0.0488 ]);
    B7	 =   exp([	 -0.5108   -0.2231 ]);
    B8	 =   exp([	 -4.6052   -0.9163 ]);
    B9   =   exp([	 -2.9957   -2.3026 ]);
    
    bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9];
    
    if(bOutputR)  
        fname = sprintf('SS_ALL_gen_%d_%i.txt',GenId,Nx);
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
end;


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
param= [1,2,3,5,6,7,8];%1:Ny-2 ;

bSynthetic = 1;

plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)