%=============================================================
%
%  Plot join plot with TMCMC samples, KDE and Marginal PDFs
%  for given data (mydata) and 2 selected parameters
%-------------------------------------------------------------
% INPUT:
%  mydata      = full matrix from curgen_db file
%  parameters  = vector of 2 parameters of interest
%  bounds      = bounds for parameters (from prior)
%  names       = names of parameters
%  groundTruth = vector of ground truth parameters
%  dump        = 0-1 to dump output
%  KDEbins     = number of bins over which to discretized kde the space
%  Nbins       = number of bins over which to compute marginal PDF
%  GenId       = Id of current generation
%  bSynthetic  = 1 for synthetic data (i.e. plot groundTruth) 0 for patient case
%=============================================================



function plotTMCMC_PairStatistic(mydata,parameters,bounds,names,groundTruth,dump,KDEbins,Nbins,GenId,bSynthetic)

% Plotting set up
%-------------------
%colors
c1 = [0.4, 0.8, 1.0]; % light blue
c2 = [0.2, 0.2, 1.0]; % dark blue
c3 = [0.6, 0.6, 1.0]; % light purple
c4 = [0.6, 0.0, 0.8]; % dark purple
c5 = [0.0, 0.8, 0.0]; % light green
c6 = [0.0, 0.2, 0.0]; % dark green
c=[c1;c2;c3;c4;c5;c6];
c=[c;c];

pointsize=10;
plotId=1;

%Data procesing parameters
N = length(parameters);

parId1=parameters(1);
parId2=parameters(2);

figure,
hold on

bestC = find( max(mydata(:,7)) == mydata(:,7));
best = mydata(bestC,:);
meanData = mean(mydata);

% I) Samples
subplot(N,N+1,plotId)
hold on
scatter(mydata(:,parId1), mydata(:,parId2),pointsize,mydata(:,end),'*','Linewidth',3);
plot(best(parId1),best(parId2),'sb','Linewidth',2);
plot(meanData(parId1),meanData(parId2),'sr','Linewidth',2);

if (bSynthetic)
    plot(groundTruth(parId1),groundTruth(parId2),'*k','Linewidth',2);
end;
xlabel( names(parId1,:) );
ylabel( names(parId2,:) );
axis( [bounds(parId1,1) bounds(parId1,2) bounds(parId2,1) bounds(parId2,2) ])
colorbar; grid on;
plotId = plotId + 1;



% II) KDE
X=linspace( bounds(parId1,1),bounds(parId1,2),Nbins);
Y=linspace( bounds(parId2,1),bounds(parId2,2),Nbins);

[JoinPDF,JFb] = getKde_2D_PDF(mydata,parameters, bounds,KDEbins,Nbins);

subplot(N,N+1,plotId)
hold on
pcolor(X,Y,JoinPDF);
if (bSynthetic)
    plot(groundTruth(parId1),groundTruth(parId2),'*k','Linewidth',2)
end;
xlabel( names(parId1,:) );
ylabel( names(parId2,:) );
 xlim(bounds(parId1,:));
ylim(bounds(parId2,:));
t=sprintf('Gaussian KDE, (ROT)');
title(t)
shading flat
colorbar
plotId = plotId + 1;

subplot(N,N+1,plotId)
hold on
pcolor(X,Y,JFb);
if (bSynthetic)
plot(groundTruth(parId1),groundTruth(parId2),'*k','Linewidth',2)
end;
xlabel( names(parId1,:) );
ylabel( names(parId2,:) );
%  xlim(bounds(parId1,:));
% ylim(bounds(parId2,:));
axis( [bounds(parId1,1) bounds(parId1,2) bounds(parId2,1) bounds(parId2,2) ])
t=sprintf('Gaussian KDE, GenId=%d',GenId);
title(t)
shading flat
colorbar
plotId = plotId + 1;



% III-IV Marginals
for parId = [parameters]
    parId
    [mPDF, bins] = getMarginalPDF(mydata,parId,bounds(parId,:),KDEbins(parId),Nbins);
    
    subplot(N,N+1,plotId)
    hold on
    z=zeros(size(mPDF));
    [fillhandle,msg]=jbfill(bins,mPDF,z,c(parId,:),c(parId,:),0,0.5);
    grid on;box on;
    xlim(bounds(parId,:));
    title(['Marginal ',names(parId,:)]);
    plotId = plotId + 1;
end

set(gcf, 'Position', [100, 100, 1049, 895]);

if(dump)
t=[names(parId1,:),'_vpcs_',names(parId2,:)];
set(gcf,'papersize',[15,15]);
print(gcf, '-djpeg', t)
end



