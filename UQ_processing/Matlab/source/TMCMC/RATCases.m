%    function PatientCases

% close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/tightfig')
addpath('../../lib/')
addpath('../../lib/freezeColors')

% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
bSynthetic = 0;
groundTruth = 1;
bOutputR = 0;
bPlotResults = 1;
bReduceUncertainties = 0;
bWriteOutput2File=0

pID = 1  % 0 synthethic

%  pathname =['../../../../RAT_DATA_F98/M',num2str(pID,'%02i'),'/Inference/'];

 % names =       ['D   ' ;'rho '   ;'ix  ';'iy  ';'iz  ';'T1uc' ;'T2uc';'Tn  '];
% groundTruth = [ 0.02, 0.6 , 0.34, 0.58, 0.48,  0.7, 0.25, 0.09  ];
% 
 pathname =['../../../../RAT_DATA_F98/M',num2str(0,'%02i'),'/Inference/'];
 names =       ['D   ' ;'rho '   ;'T1uc' ;'T2uc';'Tn  '];
% groundTruth = [ 0.02, 0.6 , 0.7, 0.25 ];


    

% pathname =['../../../../RAT_DATA_F98/M',num2str(0,'%02i'),'/Inference/'];
% names =       ['rho ' ];


if(pID == 0)
        
    bSynthetic = 1;
    groundTruth = [ 0.02, 0.6 , 0.34, 0.58, 0.48,  0.7, 0.25, 0.09  ];
    
    GenId = 22;
    myfilename = [pathname,'M00_1K_A/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    
    B0  = exp ([   -6.9078    -0.6931 ]);
    B1  = exp ([   -5.2983    0.0     ]);
    B2  = exp ([   -1.1548   -0.9957 ]);
    B3  = exp ([   -0.5946   -0.5007 ]);
    B4  = exp ([   -0.7934   -0.6799 ]);
    B5  = exp ([   -0.5108   -0.0513 ]);
    B6  = exp ([   -2.9957   -0.5108 ]);
    B7  = exp ([   -2.9957   -2.3026 ]);
    
    bounds = [B0;B1;B2;B3;B4;B5;B6;B7];
end;

if(pID == 1)
        
    bSynthetic = 1;
    groundTruth = [ 0.02, 0.6 ,0.7, 0.25, 0.05 ];
    
    GenId = 26;
    myfilename = [pathname,'M00_NEW_like_tails_5param_correc/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    
    B0  = exp ([   -6.9078    -1.3863 ]);
    B1  = exp ([   -5.2983    0.6931     ]);
    B2  = exp ([   -0.5108   -0.0513 ]);
    B3  = exp ([   -2.9957   -0.5108 ]);
    B4  = exp ([   -2.9957   -2.3026 ]);
   
    bounds = [B0;B1;B2;B3;B4];
    
end;







if(pID == 38)
            
    GenId = 21;
    myfilename = [pathname,'M38_500/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    
    B0  = exp ([   -6.9078    -1.3863 ]);
    B1  = exp ([   -5.2983    0.0     ]);
    B2  = exp ([   -1.2730   -0.6931 ]);
    B3  = exp ([   -1.4735   -0.9514 ]);
    B4  = exp ([   -1.0580   -0.6848 ]);
    B5  = exp ([   -0.5108   -0.1054 ]);
    B6  = exp ([   -2.9957   -0.5108 ]);
    B7  = exp ([   -2.9957   -2.3026 ]);
end;


if(pID == 42)
            
    GenId = 22;
    myfilename = [pathname,'M42_500/curgen_db_0',num2str(GenId),'.txt'];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = 1:Ny-1
        mydata(:,i) = exp( mydata(:,i) );
    end;
    
    
    B0  = exp ([   -6.9078    -1.3863 ]);
    B1  = exp ([   -5.2983    0.0     ]);
    B2  = exp ([   -1.1682   -0.8797 ]);
    B3  = exp ([   -1.2558   -0.9446 ]);
    B4  = exp ([   -0.9603   -0.7200 ]);
    B5  = exp ([   -0.5108   -0.0513 ]);
    B6  = exp ([   -2.9957   -0.5108 ]);
    B7  = exp ([   -2.9957   -2.3026 ]);
    
    bounds = [B0;B1;B2;B3;B4;B5;B6;B7];
end;






%% 4) Plot results - Marginals, KDE and scattered samples
%----------------------------------------------------------------


if (bPlotResults)
    
    KDEbins = 15*ones(1,Ny-1);
    Nbins   = 200;
    dump = 1;
    param= 1:Ny-1 ;
    
    plotTMCMC_AllStatistic(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)
end;





% Plot D/rho, T*rho, sigma
if(bOutputR)

        param = [1,2,3,7,Ny];
        fname = sprintf('SyntheticRhoScaled.txt');
        
    fid = fopen(fname, 'wt'); % Open for writing
    [Nx,Ny] = size(mydata);
    
    for i=1:Nx
        parm1 = mydata(i,1) / mydata(i,2); % D/rho
        parm2 = mydata(i,3) * mydata(i,2); % D*rho
        parm3 = mydata(i,7);
        parm4 = mydata(i,Ny);
        
        fprintf(fid, '%d  %d  %d  %d  \n', parm1, parm2, parm3, parm4 );

    end;
    fclose(fid);
end;

%%
if( bWriteOutput2File)
    
    outputFolder='~/software/TexmakerMacosx64/PNAS-2017-HGG/Data/';
    outputFile  = sprintf('Inference_P%d',pID);
    
    names =       ['D      ' ;'rho    '  ;'Tend   ' ;'ix     ';'iy     ';'iz     ';'PETn   ';'b      ';'T1uc   ' ;'T2uc   ';'Tn     ' ];
    
    bestC = find( max(mydata(:,end)) == mydata(:,end));
    MAPtmp = mydata(bestC,:);
    MAP = MAPtmp(1,:)
    MeanData = mean(mydata)
    varData = var(mydata)
    stdData = sqrt(varData)
    
    % write data matrix into file
    fileID = fopen([outputFolder,outputFile,'.txt'],'w');
    
    for i = 1:length(MAP)-1
        fprintf(fileID, '%4.4f  ',MAP(i));
    end;
    
    fprintf(fileID, '\n');
    
    for i = 1:length(MAP)-1
        fprintf(fileID, '%4.4f  ',MeanData(i));
    end;
    
    fprintf(fileID, '\n');
    
    for i = 1:length(MAP)-1
        fprintf(fileID, '%4.4f  ',stdData(i));
    end;
    
    fclose(fileID);
    
    
    % verbose file
    fileID = fopen([outputFolder,outputFile,'_verbose.txt'],'w');
    
    fprintf(fileID, '%s       ', 'names' );
    for i=1:length(names)
        fprintf(fileID, '%s   ', names(i,:) );
    end;
    
    fprintf(fileID, '\n');
    fprintf(fileID, '%s \n','-------------------------------------------------------------------------------------------------------------');
    fprintf(fileID, '%s    ','MAP');
    
    for i = 1:length(MAP)-1
        fprintf(fileID, '%4.4f    ',MAP(i));
    end;
    
    fprintf(fileID, '\n');
    fprintf(fileID, '%s \n','-------------------------------------------------------------------------------------------------------------');
    fprintf(fileID, '%s   ','Mean');
    
    
    for i = 1:length(MAP)-1
        fprintf(fileID, '%4.4f    ',MeanData(i));
    end;
    
    fprintf(fileID, '\n');
    fprintf(fileID, '%s \n','-------------------------------------------------------------------------------------------------------------');
    fprintf(fileID, '%s   ','std ');
    
    
    for i = 1:length(MAP)-1
        fprintf(fileID, '%4.4f    ',stdData(i));
    end;
    
    fprintf(fileID, '\n');
    fprintf(fileID, '%s \n','-------------------------------------------------------------------------------------------------------------');
    
    
    fclose(fileID);
end;




