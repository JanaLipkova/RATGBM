% function SyntheticCases

%  function SpaceSearchTests

% close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/tightfig')
addpath('../../lib/')
addpath('../../lib/freezeColors')

% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
bSynthetic      = 1;
bOutputR        = 1;
bPhysicalUnits  = 1;  % convert cm to mm
Case            = 2;  % 1 for small synthetic, 2 for big one

pathname = '../../../../../Volumes/brutus_scratch/TMCMC/';



if(Case == 0)
    
    GenId = 12;
    
    pathname = '../../../../../Volumes/brutus_scratch/';
    myfilename = sprintf('forPanos/Synthethic_CMA_4parameters/curgen_db_%03d.txt',GenId);
    
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'c   ' ;'PETn'];
    groundTruth = [ 1.3e-03 , 2.5e-02, 2280.1, 0.023 ];
    
    %     names =       ['D   ' ;'rho '  ;'T   ' ;'PETn'];
    %     groundTruth = [ 1.3e-03 , 2.5e-02, 302, 0.023 ];
    
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
%     mydata(:,3) = mydata(:,3)./mydata(:,2);
%     mydata(:,3) = sqrt(mydata(:,3));
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    %     B3  =      [     30,       1500   ];
    B3  = exp( [ 5.1416    8.7119]);
    B4  = exp ([    -4.1997   -3.2188 ]);
    
    
    if(bPhysicalUnits)
        mydata(:,1)  =   mydata(:,1) * 100;  % convert from cm^2 to mm^2
        B1 = B1 * 100;
        groundTruth(1) = groundTruth(1) * 100;
    end;
    
    bounds = [B1;B2;B3;B4];
    
end;


if(Case == 1)
    
    GenId = 11;
    %     myfilename = sprintf('SpaceSearching/SynthethicAll_All/SynthethicAll_All_4K_Lustre_C/curgen_db_%03d.txt',GenId);
    
    myfilename = sprintf('SpaceSearching/SyntheticTests/Synthetic_Movable_DrTsigma_TimeRestr/curgen_db_%03d.txt',GenId);
    
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    %     names =       ['D   ' ;'rho '  ;'Tend' ;'ix  ';'iy  ';'iz  ';'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    %     groundTruth = [ 1.3e-03 , 2.5e-02, 302 , 0.315, 0.67,   0.5,  0.019,0.8792, 0.7,   0.25,  5.0e-02  ];
    
    names =       ['D   ' ;'rho '  ;'c   ' ;'PETn'];
    groundTruth = [ 1.3e-03 , 2.5e-02, 2280.1, 0.023 ];
    
    
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1	= exp ([	-8.9480   -3.2702 ]);
    B2	= exp ([	-5.9145   -1.6607 ]);
    %     B3  =      [     30,       1500   ];
    B3  = exp( [ 5.1416    8.7119]);
    B4  = exp ([    -1.2877   -1.0261 ]);
    B5  = exp ([    -0.4677   -0.3440 ]);
    B6  = exp ([    -0.7910   -0.6238 ]);
    B7  = exp ([    -4.1997   -3.2188 ]);
    B8  = exp ([    -0.5108    0.0198 ]);
    B9  = exp ([    -0.5108   -0.1625 ]);
    B10 = exp ([    -2.9957   -0.6931 ]);
    B11 = exp ([    -2.9957   -2.3026 ]);
    
    
    if(bPhysicalUnits)
        mydata(:,1)  =   mydata(:,1) * 100;  % convert from cm^2 to mm^2
        B1 = B1 * 100;
        groundTruth(1) = groundTruth(1) * 100;
    end;
    
    bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9;B10;B11];
    
end;





if(Case == 2)
    
    GenId = 29;
    myfilename = sprintf('PatientCases/SyntheticBig/SynthethicAll_new_6K_step3/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;'rho '  ;'Tend' ;'ix  ';'iy  ';'iz  ';'PETn';'b   ';'T1uc' ;'T2uc';'Tn  '];
    groundTruth = [ 0.001715 , 0.010, 900 , 0.6, 0.45,   0.65,  0.0263, 0.7593, 0.7,   0.25,  5.0e-02  ];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1     = exp ([         -8.9480   -3.5972 ]);
    B2     = exp ([         -5.9145   -1.6607 ]);
    B3     =     ([              30    3000   ]);
    B4     = exp ([         -0.6583   -0.4023 ]);
    B5     = exp ([         -0.9789   -0.6410 ]);
    B6     = exp ([         -0.5299   -0.3016 ]);
    B7     = exp ([         -4.1997   -1.6094 ]);
    B8     = exp ([         -0.5108    0.0198 ]);
    B9     = exp ([         -0.5108   -0.2231 ]);
    B10    = exp ([         -2.9957   -0.5108 ]);
    B11    = exp ([         -2.9957   -2.3026 ]);
    
    
    if(bPhysicalUnits)
        mydata(:,1)  =   mydata(:,1) * 100;  % convert from cm^2 to mm^2
        B1 = B1 * 100;
        groundTruth(1) = groundTruth(1) * 100;
    end;
    
    bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9;B10;B11];
    
end;
%%

% 4) Plot results - Marginals, KDE and scattered samples
%----------------------------------------------------------------
KDEbins = 15*ones(1,Ny-1);
Nbins   = 200;
dump = 1;
param= 1:Ny-1 ;
% param = [1,2,3,7,8];


% rho confidence intervals
rho = mydata(:,2);
P=0.95;
confidenceIntervals(rho, P);
confidenceIntervals(rho, 0.65);


plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)

% param = [1,2,3,7,8,Ny];
param= [1,2,3,7,Ny] ;

if(bOutputR)
    fname = sprintf('SyntheticCase_%d_6K.txt',Case);
    fid = fopen(fname, 'wt'); % Open for writing
    [Nx,Ny] = size(mydata);
    for i=1:Nx
        for j=param
            fprintf(fid, '%d ', mydata(i,j));
        end
        fprintf(fid, '\n ');
    end;
    fclose(fid);
end;
