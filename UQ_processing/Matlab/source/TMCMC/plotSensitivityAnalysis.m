%=============================================================
%
%  Plot scatter plot for all samples from required generations
%
%-------------------------------------------------------------
% INPUT:
%  pathToData    = path to the data
%  GenId        = vector of generations to be plotted

%  param       = vector of parameters of interest
%  names       = names of parameters
%  bounds      = bounds for parameters (from prior)
%  groundTruth = vector of ground truth parameters
%  dump        = 0-1 to dump output
%  KDEbins     = number of bins over which to discretized kde the space
%  Nbins       = number of bins over which to compute marginal PDF
%  GenId       = Id of current generation
%  bSynthetic  = 1 for synthetic data (i.e. plot groundTruth) 0 for patient case
%=============================================================


function plotSensitivityAnalysis(pathToData, GenId, param, names, bounds, scaling )


% 1) Ploting setup
%----------------------------------------------------------------
pointsize=10;
%  b1 = [1.0, 1.0, 0.0]; % yellow
% b2 = [0.2, 0.8, 0.0]; % green
% b3 = [0.4, 0.2, 0.8]; % purple
% % b4 = [0.0, 0.0, 1.0]; % blue
% b5 = [1.0, 0.0, 0.0]; % red
% b4 = [0.0, 0.0, 0.0];  %black

b1 = [0.8, 1.0, 1.0]; % baby blue
b2 = [0.2, 0.8, 1.0]; %  light blue
b3 = [0.4, 0.2, 1.0]; % ligh purple
b4 = [1.0, 0.0, 0.0]; % red

b = [b1;b2;b3;b4];


figure,
hold on
fs=16;
set(gca,'Fontsize',fs);
sx = 6;
sy = 6;

genIndex = 1;
for nGen = GenId
    
    gen = sprintf('curgen_db_%03d.txt',nGen);
    myfilename = [pathToData,gen];
    mydata = importdata(myfilename);
    [Nx,Ny] = size(mydata);
    
    for i = param
        mydata(:,i) = exp( mydata(:,i) ) .* scaling(i);
    end;
    
    % Plot results
    %----------------------------------------------------------------
    N = length(param);
    plotId = 1;
    
    for i = 1:N
        for j = 1:N
            
            if(j > i)
                
                parId1 = param(i);
                parId2 = param(j);
                subplot(sx,sy,plotId)
                set(gca,'Fontsize',fs);
                
                hold on
                scatter(mydata(:,parId1), mydata(:,parId2),pointsize,b(genIndex,:),'*','Linewidth',3);
                
                xlabel( names(parId1,:) );
                ylabel( names(parId2,:) );
                axis( [bounds(parId1,1) bounds(parId1,2) bounds(parId2,1) bounds(parId2,2) ])
                %colorbar;
                grid on;
                
                plotId = plotId + 1;
                
            end;
        end;
    end;
    genIndex = genIndex +1;
    
end;
