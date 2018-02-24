%function SyntheticCaseParametricIC_Cartesian

close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/')

% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
bSynthetic = 1;
bOutputR   = 0;
version    = 1;  

if(version == 1)
    
    GenId      = 11;
    myfilename = sprintf('../../../TMCMC/SmallSyntheticParametricIC/SyntheticAll_Cartesian/SynthethicAll_PIC_Cart_8K/curgen_db_%03d.txt',GenId);
    mydata     = importdata(myfilename);
    [Nx,Ny]    = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'ix  ';'iy  ';'iz  '; 'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    groundTruth = [ 1.3e-03,2.5e-02, 302   , 0.315, 0.67 , 0.5  , 0.019,0.8792, 0.7,   0.25,  5.0e-02  ];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3	= exp ([	 3.4012    7.3132 ]);
    B4  = exp ([    -1.4497   -0.9172 ]);
    B5  = exp ([    -0.5358   -0.2874 ]);
    B6  = exp ([    -0.8864   -0.5496 ]);
    B7	= exp ([	-4.6052   -0.9163 ]);
    B8	= exp ([	-0.9163    0.0488 ]);
    B9	= exp ([	-0.5108   -0.2231 ]);
    B10	= exp ([	-4.6052   -0.9163 ]);
    B11	= exp ([	-2.9957   -2.3026 ]);
    
    bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9,B10;B11];
    
    if(bOutputR)
        fname = sprintf('Syn_cIC_gen_%d_%i.txt',GenId,Nx);
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





%% 4) Plot results - Marginals, KDE and scattered samples
%----------------------------------------------------------------
KDEbins = 15*ones(1,Ny-1);
Nbins   = 200;
dump = 0;

% param= 1:Ny-1 ;
param = [1,2,3,8,9,10]
plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)

param = [4,5,6]
plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)
