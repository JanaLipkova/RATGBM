%=============================================================
%
%  Plot samples from UQ inference
%-------------------------------------------------------------
%=============================================================


function PatientCases_MRIonly

addpath('../../lib/jbfill')
addpath('../../lib/tightfig')
addpath('../../lib/')
addpath('../../lib/freezeColors')
addpath('../../lib/kde2d')
addpath('../../lib/@kde')


% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
bSynthetic = 1;
bOutputR = 1;
groundTruth = 1;
bPhysicalUnits = 1;  % convert D to mm^2/day
bPlotResults = 0;


pID = 101  % 101- synt. small, 103 synt. big

pathname =['../../../../',num2str(pID,'%02i'),'/Results/Inference/'];
names =   ['D   ' ;'rho '  ;'Tend' ;'ix  ';'iy  ';'iz  ';'T1uc' ;'T2uc';'Tn  '];

if(pID == 1)
    
    GenId = 12;
    myfilename = [pathname,'P01_MRIonly_4K/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1       = exp ([     -8.9480   -3.2702  ]);
    B2       = exp ([     -5.9145   -1.6607  ]);
    %   B3       = exp ([      5.1416    8.9623  ]);  % c
    B3       =      [          30,    1700    ];   % Tend
    B4       = exp ([     -0.4992   -0.3515  ]);
    B5       = exp ([     -0.4513   -0.3101  ]);
    B6       = exp ([     -0.6919   -0.5154  ]);
    B7       = exp ([     -0.5108   -0.2231  ]);
    B8       = exp ([     -2.9957   -0.5108  ]);
    B9       = exp ([     -2.9957   -2.3026  ]);
end;


if(pID == 6)
    
    GenId = 18;
    myfilename = sprintf('Patient06/P06_step3_6K/curgen_db_%03d.txt',GenId);
    
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1     = exp ([     -9.2103   -3.2702  ]);
    B2     = exp ([     -6.2146   -1.6607  ]);
    B3     =    ([           60,   2900    ]);
    B4      = exp ([     -1.0921   -0.5789 ]);
    B5      = exp ([     -0.6587   -0.2977 ]);
    B6      = exp ([     -0.6387   -0.3175 ]);
    B7     = exp ([     -4.0174   -1.2040  ]);
    B8     = exp ([     -0.5108    0.0198  ]);
    B9     = exp ([     -0.6931   -0.2231  ]);
    B10    = exp ([     -2.9957   -0.5108  ]);
    B11    = exp ([     -2.9957   -2.3026  ]);
end;


if(pID == 7)
    
    GenId = 13;
    myfilename = [pathname,'P07_MRIonly_4K/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1     = exp ([     -9.2103   -3.2702  ]);
    B2     = exp ([     -6.2146   -1.6607  ]);
    B3     =     ([           30,   1800   ]);
    B4     = exp ([     -0.9903   -0.7476  ]);
    B5     = exp ([     -0.5592   -0.3950  ]);
    B6     = exp ([     -0.4420   -0.2946  ]);
    B7     = exp ([     -0.6931   -0.2231  ]);
    B8     = exp ([     -2.9957   -0.5108  ]);
    B9     = exp ([     -2.9957   -2.3026  ]);
end;


if(pID == 9)
    
    GenId = 12;
    
    myfilename = [pathname,'P09_MRIonly_4K/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1   = exp ([   -8.9480   -3.2702  ]);
    B2   = exp ([   -5.9145   -1.6607  ]);
    B3   =     ([        30,   1400    ]);
    B4   = exp ([   -1.4271   -1.0788  ]);
    B5   = exp ([   -0.5025   -0.3216  ]);
    B6   = exp ([   -0.6733   -0.4943  ]);
    B7   = exp ([   -0.6931   -0.2231  ]);
    B8  = exp ([   -2.9957   -0.5108  ]);
    B9  = exp ([   -2.9957   -2.3026  ]);
end;


if(pID == 11)
    
    GenId = 19;
    myfilename = sprintf('Patient11/P11_step3_6K/curgen_db_%03d.txt',GenId);
    
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1     = exp ([        -8.9480   -3.2702 ]);
    B2     = exp ([        -5.9145   -1.6607 ]);
    B3     =     ([            60     3100   ]);
    B4     = exp ([        -1.2730   -0.7765 ]);
    B5     = exp ([        -0.8153   -0.3603 ]);
    B6     = exp ([        -0.4620   -0.2107 ]);
    B7     = exp ([        -4.0174   -1.5141 ]);
    B8     = exp ([        -0.5108    0.0198 ]);
    B9     = exp ([        -0.5108   -0.2231 ]);
    B10    = exp ([        -2.9957   -0.5108 ]);
    B11    = exp ([        -2.9957   -2.3026 ]);
end;

if(pID == 12)
    
    GenId = 16;
    myfilename = sprintf('Patient12/P12_step3_6K/curgen_db_%03d.txt',GenId);
    
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1     = exp ([          -8.9480   -3.2702 ]);
    B2     = exp ([         -5.9145   -1.6607 ]);
    % B3     = exp ([           6.5280   10.1638 ]);
    B3     =     ([             60     2900  ]);
    B4     = exp ([        -0.8675   -0.5978 ]);
    B5     = exp ([        -0.8916   -0.4943 ]);
    B6     = exp ([        -0.3496   -0.1567 ]);
    B7     = exp ([        -4.1997   -1.3863 ]);
    B8     = exp ([        -0.5108    0.0198 ]);
    B9     = exp ([        -0.5108   -0.2231 ]);
    B10    = exp ([        -2.9957   -0.5108 ]);
    B11    = exp ([        -2.9957   -2.3026 ]);
end;



if(pID == 20)
    
    GenId = 19;
    myfilename = sprintf('Patient20/P20_NC_step3_6K_sigma_025_wait1500_newAnatomy/curgen_db_%03d.txt',GenId);
    
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    
    B1     = exp ([          -8.9480   -3.2702 ]);
    B2     = exp ([         -5.9145   -1.6607 ]);
    % B3     = exp ([           6.5280   10.1638 ]);
    B3     =     ([              60     3100  ]);
    B4     = exp ([         -1.1441   -0.6568 ]);
    B5     = exp ([         -1.1317   -0.6116 ]);
    B6     = exp ([         -0.7731   -0.4595 ]);
    B7     = exp ([        -4.0174   -1.3863 ]);
    B8     = exp ([        -0.5108    0.0198 ]);
    B9     = exp ([        -0.5108   -0.2231 ]);
    B10    = exp ([         -2.9957   -0.5108 ]);
    B11    = exp ([        -2.9957   -2.3026 ]);
end;


if(pID == 22)
    
    GenId = 18;
    myfilename = sprintf('Patient22/V2/P22_V2_new_6K_sigma_022/curgen_db_%03d.txt',GenId);
    
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1     = exp ([   -8.9480   -3.2702  ]);
    B2     = exp ([   -5.9145   -1.6607  ]);
    %     B3     = exp ([    5.1417   10.0982  ]);  %c
    B3     =     ([        30, 3000      ]);
    B4     = exp ([   -0.5874   -0.3563  ]);
    B5     = exp ([   -0.9127   -0.6052  ]);
    B6     = exp ([   -0.5430   -0.3208  ]);
    B7     = exp ([   -4.0174   -1.3863  ]);
    B8     = exp ([   -0.5108    0.0198  ]);
    B9     = exp ([   -0.6931   -0.2231  ]);
    B10    = exp ([   -2.9957   -0.5108  ]);
    B11    = exp ([   -2.9957   -2.3026  ]);
end;


if(pID == 23)
    
    GenId = 15;
    myfilename = [pathname,'Patient23_MRI_4K/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    
    B1   = exp ([           -8.9480   -3.2702 ]);
    B2   = exp ([           -5.9145   -1.6607 ]);
    %     B3   = exp ([            6.5280    9.8875 ]);
    B3   =     ([                60,    3900   ]);
    B4   = exp ([            -0.9808   -0.5025 ]);
    B5   = exp ([           -0.5621   -0.2744  ]);
    B6   = exp ([           -0.5709   -0.3079  ]);
    B7   = exp ([           -0.5108   -0.2231  ]);
    B8   = exp ([           -2.9957   -0.5108  ]);
    B9   = exp ([           -2.9957   -2.3026  ]);
end;



if(pID == 31)
    
    GenId = 23;
    
    myfilename = [pathname,'Patient31_MRI_4K/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    
    B1     = exp ([   -8.9480   -3.2702  ]);
    B2     = exp ([   -5.9145   -1.6607  ]);
    B3     =     ([        60,   2700    ]);
    B4     = exp ([   -0.9808   -0.5025 ]);
    B5     = exp ([   -0.5621   -0.2744 ]);
    B6     = exp ([   -0.5709   -0.3079 ]);
    B7     = exp ([   -0.5108   -0.2231 ]);
    B8     = exp ([   -2.9957   -0.5108 ]);
    B9     = exp ([   -2.9957   -2.3026 ]);
    
    
end;


% Synthetic cases
if(pID == 101)
    
    bSynthetic = 1;
    groundTruth = [ 1.3e-03 , 2.5e-02, 302 , 0.315, 0.67,   0.5, 0.7,   0.25,  5.0e-02  ];
    
    GenId = 17;
    pathname ='../../../../S1/Results/Inference/';
    myfilename = [pathname,'Synthetic_MRIonly_4K/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1   = exp ([    -8.9480   -3.2702 ]);
    B2   = exp ([    -5.9145   -1.6607 ]);
    %     B3   = exp ([     5.1416    8.7119 ]);
    B3       =  [          30,    1500 ];
    B4   = exp ([    -1.2877   -1.0261 ]);
    B5   = exp ([    -0.4677   -0.3440 ]);
    B6   = exp ([    -0.7910   -0.6238 ]);
    B7   = exp ([    -0.5108   -0.1625 ]);
    B8   = exp ([    -2.9957   -0.6931 ]);
    B9   = exp ([    -2.9957   -2.3026 ]);
    
end;

% Synthetic case
if(pID == 103)
    
    bSynthetic = 1;
    groundTruth = [ 0.001715 , 0.010, 900 , 0.6, 0.45,   0.65, 0.8,   0.3,  5.0e-02  ];
    
    
    GenId = 18;
    pathname ='../../../../S103/Results/Inference/';
    myfilename = [pathname,'SyntheticBig_P103_MRI_4K_newPrior/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    mydata(:,3) = mydata(:,3)./mydata(:,2);
    mydata(:,3) = sqrt(mydata(:,3));
    
    B1      = exp ([   -8.9480   -3.5972 ]);
    B2      = exp ([   -5.9145   -1.6607 ]);
    %     B3      = exp ([    5.1417   10.0982 ]);
    B3      =  [          30,    3000 ];
    B4      = exp ([   -0.6583   -0.4023 ]);
    B5      = exp ([   -0.9789   -0.6410 ]);
    B6      = exp ([   -0.5299   -0.3016 ]);
    B7      = exp ([   -0.5108   -0.1054 ]);
    B8      = exp ([   -2.9957   -0.5108 ]);
    B9      = exp ([   -2.9957   -2.3026 ]);
    
end;


%% Physical units

if(bPhysicalUnits)
    mydata(:,1)  =   mydata(:,1) * 100;  % convert from cm^2 to mm^2
    B1 = B1 * 100;
    
    if(bSynthetic)
        groundTruth(1) = groundTruth(1).*100;
    end;
end;

bounds = [B1;B2;B3;B4;B5;B6;B7;B8;B9];

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
confidenceIntervals(rho, 0.95);
confidenceIntervals(rho, 0.65);

D = mydata(:,1);
alpha = D./rho;
mean_aplha = mean(alpha(:))
std_alpha = std(alpha(:))

if (bPlotResults)
    plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)
end;




if(bOutputR)
    
    if(bSynthetic)
        param = [1,2,3,Ny];
        fname = sprintf('SynthethicSmallMRIonly_4K.txt');
    else
        param = [1,2,3,7,Ny];
        fname = sprintf('Patient_%d_6K_final.txt',Patient);
    end;
    
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
