% function SpineSyntheticCases

addpath('../../lib/jbfill')
addpath('../../lib/tightfig')
addpath('../../lib/')

% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
bSynthetic      = 1;
bOutputR        = 0;
bPhysicalUnits  = 0;  % convert cm to mm
Case            = 5;  % 1 for small synthetic, 2 for big one

pathname = '../../../../../Volumes/brutus_scratch/TMCMC/SpineTumours/Synthetic/';


if(Case == 1)
    
    GenId = 26;
    myfilename = sprintf('Model1/M1_v23_2K/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;  'rho '  ;'Tend' ;'c   ';'ix  ';'iy  ';'iz  ';'PETn';'b   '];
    groundTruth = [ 1.3e-03, 1.1e-02,  730,     0.1,   0.32,   0.42, 0.65,  0.0260, 0.7687];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    B0   = exp ([  -10.5197   -5.8091]);
    B1   = exp ([  -8.2171   -3.4420 ]);
    B2   = exp ([   3.4012    7.3132 ]);
    B3   = exp ([  -4.6052    0.0198 ]);
    B4   = exp ([  -1.5572   -0.9145 ]);
    B5   = exp ([  -1.1251   -0.6643 ]);
    B6   = exp ([  -0.6138   -0.3129 ]);
    B7   = exp ([  -4.1997   -2.3026 ]);
    B8   = exp ([  -1.2040    0.0198 ]);
    
    bounds = [B0;B1;B2;B3;B4;B5;B6;B7;B8];
end;


if(Case == 2)
    
    GenId = 31;
    myfilename = sprintf('Model1/M1_v23_fixedIC_2K/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['D   ' ;  'rho '  ;'Tend' ;'c   ';'PETn';'b   '];
    groundTruth = [ 1.3e-03, 1.1e-02,  730,   0.1,   0.0260, 0.7687];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    B0   = exp ([  -10.5197   -5.8091]);
    B1   = exp ([  -8.2171   -3.4420 ]);
    B2   = exp ([   3.4012    7.3132 ]);
    B3   = exp ([  -4.6052    0.0198 ]);
    B4   = exp ([  -4.1997   -2.3026 ]);
    B5   = exp ([  -1.2040    0.0198 ]);
    
    bounds = [B0;B1;B2;B3;B4;B5];
    
end;


if(Case == 3)
    
    GenId = 34;
    myfilename = sprintf('Model2/M2_v23_2K/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['d   ' ;  'rho '  ;'Tend';'ix  ';'iy  ';'iz  '; 'PETn';'b   '];
    groundTruth = [ 0.247, 1.1e-02,   730,   0.32, 0.42, 0.65, 0.03, 0.6336];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    B0  = exp ([  -5.9522    1.2809 ]);
    B1  = exp ([  -8.2171   -3.4420 ]);
    B2  = exp ([   3.4012    7.3132 ]);
    B3  = exp ([  -1.4679   -0.8906 ]);
    B4  = exp ([  -1.0939   -0.6638 ]);
    B5  = exp ([  -0.6037   -0.3191 ]);
    B6  = exp ([  -4.1997   -2.3026 ]);
    B7  = exp ([  -1.2040    0.0198 ]);
    
    bounds = [B0;B1;B2;B3;B4;B5;B6;B7];
end;

if(Case == 4)
    
    GenId = 21;
    myfilename = sprintf('Model2/M2_v23_fixedIC_2K/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['d   ' ;  'rho '  ;'Tend';'PETn';'b   '];
    groundTruth = [ 0.247, 1.1e-02,    730,  0.03, 0.6336];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
B0  = exp ([  -5.9522    1.2809 ]);
B1  = exp ([  -8.2171   -3.4420 ]);
B2  = exp ([   3.4012    7.3132 ]);
B3  = exp ([  -4.1997   -2.3026 ]);
B4  = exp ([  -1.2040    0.0198 ]);
    
    bounds = [B0;B1;B2;B3;B4];
end;


if(Case == 5)
    
    GenId = 15;
    myfilename = sprintf('Model2/Age_70/M2_70_v23_1K/curgen_db_%03d.txt',GenId);
    myfilename = [pathname,myfilename];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    names =       ['d   ' ;  'rho '  ;'Tend';'ix  ';'iy  ';'iz  '; 'PETn';'b   '];
    groundTruth = [ 0.08, 0.014,   730,   0.41, 0.265, 0.32, 0.0245, 0.76854];
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    B0   = exp ([  -4.8036    0.9933 ]);
    B1   = exp ([  -7.1039   -3.6119 ]);
    B2   = exp ([   3.4012    7.3132 ]);
    B3   = exp ([  -1.2275   -0.7486 ]);
    B4   = exp ([  -1.6591   -0.9934 ]);
    B5   = exp ([  -1.4500   -0.8805 ]);
    B6   = exp ([  -4.1997   -2.3026 ]);
    B7   = exp ([  -0.9163    0.0198 ]);
    
    bounds = [B0;B1;B2;B3;B4;B5;B6;B7];
end;

%%

% 4) Plot results - Marginals, KDE and scattered samples
%----------------------------------------------------------------
KDEbins = 15*ones(1,Ny-1);
Nbins   = 200;
dump = 1;
param= 1:Ny-1 ;


plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)

% param = [1,2,3,7,8,Ny];
param= [1,2,3,7,Ny] ;

if(bOutputR)
    fname = sprintf('SyntheticCase_%d_4K.txt',Case);
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
