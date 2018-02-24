% function ParametricCorrelationsTest

% close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/tightfig')
addpath('../../lib/')
addpath('../../lib/freezeColors')

% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
bSynthetic = 1;
bOutputR = 0;
version = 2;  % 1) rho = m*D, 2) rho = c/T^2, 3) rho = m*D = c / T^2

if(version == 1)
    
    GenId = 14;
    myfilename = sprintf('../../../TMCMC/SmallSynthetic/SynthethicAll_manifolds/Synthetic_rho_mD/SynthethicAll_mD_12K/curgen_db_%03d.txt',GenId);
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
     
    names =       ['D   ' ;'rho '  ;'Tend' ;'ICn '; 'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302   , 0.0062, 0.019,0.8792, 0.7,   0.25,  5.0e-02  ];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,1) = mydata(:,2)./mydata(:,1);
    
    B1 =       [    1.3e-04,    3.8e-02];
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3	= exp ([	 3.4012    7.3132 ]);
    B4	= exp ([	-6.9078   -1.2040 ]);
    B5	= exp ([	-4.6052   -0.9163 ]);
    B6	= exp ([	-0.9163    0.0488 ]);
    B7	= exp ([	-0.5108   -0.2231 ]);
    B8	= exp ([	-4.6052   -0.9163 ]);
    B9	= exp ([	-2.9957   -2.3026 ]);
    
    bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9];
    
    if(bOutputR)
        fname = sprintf('Syn_mD_gen_%d_%i.txt',GenId,Nx);
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
if(version == 2)
    
    GenId = 13;
    
    
    myfilename = sprintf('../../../TMCMC/SmallSynthetic/SynthethicAll_manifolds/Synthetic_rho_ciT_reduced/SynthethicAll_ciT_r_10K/curgen_db_%03d.txt',GenId);
    %     myfilename = sprintf('../../../TMCMC/SmallSynthetic/SynthethicAll_manifolds/Synthetic_rho_ciT/SynthethicAll_ciT_12K/curgen_db_%03d.txt',GenId);
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);

    names =       ['D   ' ;'rho '  ;'Tend' ;'ICn '; 'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302   , 0.0062, 0.019,0.8792, 0.7,   0.25,  5.0e-02  ];
        
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
     mydata(:,3) = mydata(:,3)./mydata(:,2);
     mydata(:,3) = sqrt(mydata(:,3));
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  =      [     30,       1500   ];
    B4	= exp ([	-6.9078   -1.2040 ]);
    B5	= exp ([	-4.6052   -0.9163 ]);
    B6	= exp ([	-0.9163    0.0488 ]);
    B7	= exp ([	-0.5108   -0.2231 ]);
    B8	= exp ([	-4.6052   -0.9163 ]);
    B9	= exp ([	-2.9957   -2.3026 ]);
    
    bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9];
    
    if(bOutputR)
        fname = sprintf('Syn_ciT2_gen_%d_%i.txt',GenId,Nx);
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
if(version == 3)
    
    GenId = 13;
    myfilename = sprintf('../../../TMCMC/SmallSynthetic/SynthethicAll_manifolds/Synthetic_rho_ciT_mD/Synthetic_rho_ciT_mD_8K/curgen_db_%03d.txt',GenId);
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
%     names =       ['m   ' ;'rho '  ;'c   ' ;'ICn '; 'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
%     groundTruth = [19.23, 2.5e-02,  2280.1 , 0.0062, 0.019,0.8792, 0.7,   0.25,  0.01  ];

    names =       ['D   ' ;'rho '  ;'Tend' ;'ICn '; 'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302   , 0.0062, 0.019,0.8792, 0.7,   0.25,  5.0e-02  ];

    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,1) = mydata(:,2)./mydata(:,1);
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1 =       [    1.3e-04,    3.8e-02];
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  =      [     30,       1500   ];
    B4	= exp ([	-6.9078   -1.2040 ]);
    B5	= exp ([	-4.6052   -0.9163 ]);
    B6	= exp ([	-0.9163    0.0488 ]);
    B7	= exp ([	-0.5108   -0.2231 ]);
    B8	= exp ([	-4.6052   -0.9163 ]);
    B9	= exp ([	-2.9957   -2.3026 ]);
    
    bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9];
    
    if(bOutputR)
        fname = sprintf('Syn_mc_gen_%d_%i.txt',GenId,Nx);
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
param= 1:Ny-1 ;
%param = [1,2,3,6,7,8]
 %param = [1,2,3]

plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)
