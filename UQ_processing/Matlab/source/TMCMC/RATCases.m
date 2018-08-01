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

pID = 38  % 0 synthethic

 pathname =['../../../../RAT_DATA_F98/M',num2str(pID,'%02i')];
 names =       ['D   ' ;'rho ' ;'s   ';'T1uc' ;'T2uc';'Tn  '];


if(pID == 0)
        
    bSynthetic = 1;
    groundTruth = [ 0.02, 0.6 , 0.806329,  0.7, 0.25, 0.09  ];
    
    GenId = 19;
    myfilename = [pathname,'/Inference/M00_500_2timeLike/curgen_db_0',num2str(GenId),'.txt'];
end;


if(pID == 1)
        
    bSynthetic = 1;
    groundTruth = [ 0.015, 0.5 , 0.795171,  0.75, 0.25, 0.05  ];
    GenId = 18;
    myfilename = [pathname,'/Inference/M01_3K_2TP/curgen_db_0',num2str(GenId),'.txt'];
end;

if(pID == 2)
        
    bSynthetic = 1;
    groundTruth = [ 0.02, 0.6 , 0.821514,  0.7, 0.25, 0.05  ];
    
    GenId = 18;
    myfilename = [pathname,'/Inference/M02_3K_2TP/curgen_db_0',num2str(GenId),'.txt'];
end;


if(pID == 29)
    GenId = 22;
    myfilename = [pathname,'/Results/Inference/M29_500_2TP/curgen_db_0',num2str(GenId),'.txt'];
end;

if(pID == 30)
    GenId = 21;
    myfilename = [pathname,'/Results/Inference/M30_2K_2TP/curgen_db_0',num2str(GenId),'.txt'];
end;

if(pID == 34)
    GenId = 23;
    myfilename = [pathname,'/Results/Inference/M34_500_2TP/curgen_db_0',num2str(GenId),'.txt'];
end;


if(pID == 36)
    GenId = 20;
    myfilename = [pathname,'/Results/Inference/M36_500_2TP/curgen_db_0',num2str(GenId),'.txt'];
end;

if(pID == 38)
    GenId = 20;
    myfilename = [pathname,'/Results/Inference/M38_500_2TP/curgen_db_0',num2str(GenId),'.txt'];
end;

if(pID == 42)
    GenId = 26;
    myfilename = [pathname,'/Results/Inference/M42_500_2TP/curgen_db_0',num2str(GenId),'.txt'];
end;

mydata = importdata(myfilename);
[Nx,Ny] = size(mydata);

for i = 1:Ny-1
    mydata(:,i) = exp( mydata(:,i) );
end;

B0   = exp ([   -6.9078    0.9163 ]);
B1   = exp ([   -5.2983    0.6931 ]);
B2   = exp ([   -0.5108    0.0    ]);
B3   = exp ([   -0.5108   -0.2231 ]);
B4   = exp ([   -2.9957   -0.5108 ]);
B5   = exp ([   -2.9957   -2.3026 ]);
 
bounds = [B0;B1;B2;B3;B4;B5];

%% 4) Plot results - Marginals, KDE and scattered samples
%----------------------------------------------------------------


if (bPlotResults)
    
%     KDEbins = 15*ones(1,Ny-1);
    KDEbins = 20*ones(1,Ny-1);
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




