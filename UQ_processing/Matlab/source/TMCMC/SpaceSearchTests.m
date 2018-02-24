%  function SpaceSearchTests

% close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/tightfig')
addpath('../../lib/')
addpath('../../lib/freezeColors')

% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
bSynthetic = 1;
bOutputR = 0;
version = 7;  % 0) D,r,C 1) D,r,T, 2) D,r,c,IC 3) D,r,c,b 4) D,r,c,b,PETsigma 5)D,r,c,T1uc,T2uc,Tis2 or same but Tis2 fixed
% 6) all moveable expect IC 7) all moveable 8) PET only with fixed IC   9) PET only all moveable

pathname = '../../../../../Volumes/brutus_scratch/TMCMC/SpaceSearching/';

if(version == 1)
    
    GenId = 11;
    myfilename = sprintf('/SynthethicAll_DcT/SynthethicAll_ciT_1K_UsingLuster_b/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];
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
    
end;

if(version == 1)
    
    GenId = 12;
    myfilename = sprintf('SynthethicAll_DrT/SynthethicAll_DrT_1K_UsingLuster_c/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'ICn '; 'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302   , 0.0062, 0.019,0.8792, 0.7,   0.25,  5.0e-02  ];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  = exp ([     3.4012    7.3132 ]);
    
    
    bounds = [B1;B2;B3];
    
end;

if(version == 2)
    
    GenId = 16;
    myfilename = sprintf('SynthethicAll_DcT_cIC/SynthethicAll_PIC_Cart_1K_Lustre/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'ix  '; 'iy  ';'iz  '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302 , 0.315, 0.67,   0.5 ];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  =      [     30,       1500   ];
    B4  = exp ([    -1.2877   -1.0261 ]);
    B5  = exp ([    -0.4677   -0.3440 ]);
    B6  = exp ([    -0.7910   -0.6238 ]);
    
    bounds = [B1;B2;B3;B4;B5;B6];
    
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


if(version == 3)
    GenId = 12;
    myfilename = sprintf('../../../TMCMC/SpaceSearching/SynthethicAll_DcT_b/SynthethicAll_ciT_b_2K/curgen_db_%03d.txt',GenId);
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'b   '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302 , 0.85 ];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  =      [     30,       1500   ];
    B4  = exp ([    -0.9163    0.0488 ]);
    
    bounds = [B1;B2;B3;B4];
    
end;


if(version == 4)
    GenId = 13;
    myfilename = sprintf('../../../TMCMC/SpaceSearching/SynthethicAll_DcT_b_PETnoise/SynthethicAll_ciT_b_PETnoise_4K/curgen_db_%03d.txt',GenId);
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend'; 'PETn';'b   '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302 , 0.019,0.8792];
    
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  =      [     30,       1500   ];
    B4   = exp ([   -4.6052   -0.9163 ]);
    B5   = exp ([   -0.9163    0.0488 ]);
    bounds = [B1;B2;B3;B4;B5];
    
end;


if(version == 5)
    GenId = 18;
    
    myfilename = sprintf('/SynthethicAll_DcT_Ti_Tsigma/SynthethicAll_ciT_Ti_Tsigma_1K_Lustre/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];    
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'T1uc' ;'T2uc';'Tn  '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302 , 0.7,   0.25,  5.0e-02  ];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  =      [     30,       1500   ];
    B4	= exp ([	-0.5108   -0.2231 ]);
    B5	= exp ([	-4.6052   -0.9163 ]);
    B6	= exp ([	-2.9957   -2.3026 ]);
    
    bounds = [B1;B2;B3;B4;B5;B6];
    
end;



if(version == 6)
    
    GenId = 20;
    myfilename = sprintf('SynthethicAll_DcT_FixIC/SynthethicAll_FixIC_1K_Lustre/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];    
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ; 'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302 , 0.019,0.8792, 0.7,   0.25,  5.0e-02  ];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  =      [     30,       1500   ];
    B4  = exp ([    -4.1997   -3.2188 ]);
    B5  = exp ([    -0.5108    0.0198 ]);
    B6  = exp ([    -0.5108   -0.2231 ]);
    B7  = exp ([    -2.9957   -0.9162 ]);
    B8  = exp ([    -2.9957   -2.3026 ]);
    
    bounds = [B1;B2;B3;B4;B5;B6;B7;B8];
end;

if(version == 7)
    
    GenId = 22;
    myfilename = sprintf('SynthethicAll_All/SynthethicAll_All_4K_Lustre_C/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'ix  ';'iy  ';'iz  ';'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302 , 0.315, 0.67,   0.5,  0.019,0.8792, 0.7,   0.25,  5.0e-02  ];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
%     B1	= exp ([	-8.9480   -3.2702 ]);
%     B2	= exp ([	-5.9145   -1.6607 ]);
%     B3  =      [     30,       1500   ];
%     B4  = exp ([    -1.2877   -1.0261 ]);
%     B5  = exp ([    -0.4677   -0.3440 ]);
%     B6  = exp ([    -0.7910   -0.6238 ]);
%     B7  = exp ([    -4.1997   -3.2188 ]);
%     B8  = exp ([    -0.5108    0.0198 ]);
%     B9  = exp ([    -0.5108   -0.0513 ]);
%     B10 = exp ([    -2.9957   -0.5108 ]);
%     B11 = exp ([    -2.9957   -2.3026 ]);
    
    

    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  =      [     30,       1500   ];
    B4  = exp ([    -1.2877   -1.0261 ]);
    B5  = exp ([    -0.4677   -0.3440 ]);
    B6  = exp ([    -0.7910   -0.6238 ]);
    B7  = exp ([    -4.1997   -1.6094 ]);
    B8  = exp ([    -0.5108    0.0198 ]);
    B9  = exp ([    -0.5108   -0.1625 ]);
    B10 = exp ([    -2.9957   -0.6931 ]);
    B11 = exp ([    -2.9957   -2.3026 ]);


    
    
    bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9;B10;B11];
    
end;

if(version == 8)
    
    GenId = 12;
    myfilename = sprintf('../../../TMCMC/SpaceSearching/SynthethicPET/SyntheticPET_FixIC/SynthethicPET_FixIC_8K/curgen_db_%03d.txt',GenId);
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'PETn';'b   ' ];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302 ,0.019,0.8792  ];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  =      [     30,       1500   ];
    B4  = exp ([    -4.6052   -0.9163 ]);
    B5  = exp ([    -0.9163    0.0488 ]);
    
    bounds = [B1;B2;B3;B4;B5];
    
end;


if(version == 9)
    
    GenId = 16;
    myfilename = sprintf('../../../TMCMC/SpaceSearching/SynthethicPET/SyntheticPET_All/SynthethicPET_All_10K/curgen_db_%03d.txt',GenId);
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'ix  ';'iy  ';'iz  ';'PETn';'b   '];
    groundTruth = [ 1.3e-03 , 2.5e-02, 302 , 0.315, 0.67,  0.5,  0.019, 0.8792];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    B3  =      [     30,       1500   ];
    B4  = exp ([    -1.2877   -1.0261 ]);
    B5  = exp ([    -0.4677   -0.3440 ]);
    B6  = exp ([    -0.7910   -0.6238 ]);
    B7  = exp ([    -4.1997   -3.2188 ]);
    B8  = exp ([    -0.5108    0.0198 ]);
    
    bounds = [B1;B2;B3;B4;B5;B6;B7;B8];
    
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
% param = [5,6]

plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)
