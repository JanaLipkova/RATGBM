% function PETvsMeanStd

addpath('../../lib/jbfill')
addpath('../../lib/tightfig')
addpath('../../lib/')
addpath('../../lib/vi')
addpath('../../lib/freezeColors')


% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
pathname ='../../../';

fs = 16;
Patient = 1;
pID = '01';
myfilename = [pathname,pID,'/Results/ResultsNii/'];

InPET  = [myfilename,'P',pID,'_PET_out.nii'];
InT1s  = [myfilename,'P',pID,'_T1wsegm_out.nii'];
InT2s  = [myfilename,'P',pID,'_FLAIRsegm_out.nii'];
InMean = [myfilename,'P',pID,'_Mean_out.nii'];
InStd  = [myfilename,'P',pID,'_Std_out.nii'];
InMAP  = [myfilename,'P',pID,'_MAP_out.nii'];


% read in data
out = MRIread(InPET);
PET = out.vol;
out = MRIread(InT1s);
T1s = out.vol;
out = MRIread(InT2s);
T2s = out.vol;
out = MRIread(InMean);
Mean = out.vol;
out = MRIread(InStd);
Std = out.vol;
out = MRIread(InMAP);
MAP = out.vol;

%%
% color maps
a1 = [0.0, 0.2, 0.0]; % dark green
b1 = [0.0, 0.8, 0.0]; % light green

a2 = [0.0, 0.2, 1.0];  % dark blue
b2 = [0.2, 0.8, 1.0];  % light blue

a3 = [0.0, 0.8, 1.0];  % dark cyan
b3 = [0.4, 1.0, 1.0];  % light cyan

a4 = [0.6, 0.0, 0.8]; % dark purple
b4 = [0.6, 0.6, 1.0]; % light purple

a5 = [1.0, 0.0, 0.8];  % dark pink
b5 = [1.0, 0.6, 1.0];  % light pink

a6 = [1.0, 0.4, 0.2]; % orange
b6 = [1.0, 1.0, 0.4]; % yellow






% Patient 01
if(Patient == 1)
    sliceZ = 140;
    sliceX = 169;
    bmean = 0.95;%12963188920;
    
    % 1) load volume data: PET, tumour segm, mean + std
    myfilename = [pathname,pID,'/Results/ResultsNii/'];
    
    
    
    
    PET01 = loadMatrix(myfilename);
%     PET01 = PET01.*bmean;
    
    myfilename = [pathname,'01/Results/Propagations/PropagationsP01_4K/PropagationMean_4319.dat'];
    Mean01 = loadMatrix(myfilename);
    
    myfilename = [pathname,'01/Results/Propagations/PropagationsP01_4K/PropagationVar_4319.dat'];
    Var01 = loadMatrix(myfilename);
    Std01 = sqrt(Var01);
    clear Var01;
    
    UB1  = Mean01(sliceX,:,sliceZ) + 2 * Std01(sliceX,:,sliceZ);
    LB1  = Mean01(sliceX,:,sliceZ) - 2 * Std01(sliceX,:,sliceZ);
    PET1 = PET01(sliceX,:,sliceZ);
         PET1 = PET1./max(PET01(:));
         PET1 = PET1.*bmean;
    Mean1s = Mean01(sliceX,:,sliceZ);
    
    x=1:256;
    
    %     figure;
    %     hold on
    %     set(gca,'Fontsize',fs);
    %     title('Patient01')
    %     plot(x,UB1')
    %     plot(x,LB1')
    %     plot(x,Mean1s')
    %     plot(x,PET1')
    %     legend('UB','LB','Mean','PETmean','Location','northwest')
    %     axis([120,200,0,1])
    %     box on
    %     grid on
    
    figure;
    hold on;
    title('Patient01')
    set(gca,'Fontsize',fs);
    jbfill(x,LB1,UB1,b2,a2,0,0.5);
    plot(x,PET1','Linewidth',2)
    plot(x,Mean1s','Linewidth',2)
    legend('2std','PET','Mean','Location','northwest')
    axis([130,210,0,1])
    box on
    grid on
    print('P01','-djpeg')
    print('P01','-dpdf')
end
%% Patient07
Patient = 7;
if(Patient == 7)
    sliceZ = 179;
    sliceY = 156;
    sliceX = 105;
    bmean =0.85;
    
    % 1) load volume data: PET, tumour segm, mean + std
    myfilename = [pathname,'07/Results/ProcessingDataHR/PET_data.dat'];
    PET07 = loadMatrix(myfilename);
    %      PET07 = PET07.*bmean;
    
    myfilename = [pathname,'07/Results/Propagations/PropagationsP07_4K/PropagationMean_4319.dat'];
    %     myfilename = [pathname,'07/Results/Propagations/PropagationsP07_V2_All_4K_old/PropagationMean_4320.dat'];
    Mean07 = loadMatrix(myfilename);
    
    %     myfilename = [pathname,'07/Results/Propagations/PropagationsP07_V2_All_4K_old/PropagationVar_4320.dat'];
    myfilename = [pathname,'07/Results/Propagations/PropagationsP07_4K/PropagationVar_4319.dat'];
    Var07 = loadMatrix(myfilename);
    Std07 = sqrt(Var07);
    clear Var07;
    
    UB7  = Mean07(:,sliceY,sliceZ) + 2 * Std07(:,sliceY,sliceZ);
    LB7  = Mean07(:,sliceY,sliceZ) - 2 * Std07(:,sliceY,sliceZ);
    PET7 = PET07(:, sliceY,sliceZ);
    PET7 = PET7./max(PET7);
    PET7 = PET7 * bmean;
    Mean7s = Mean07(:,sliceY,sliceZ);
    
    %%
    
    x=1:256;
    
    %     figure;
    %     hold on
    %     set(gca,'Fontsize',fs);
    %     title('Patient07')
    %     plot(x,UB7','Linewidth',2)
    %      plot(x,LB7','Linewidth',2)
    %     plot(x,Mean7s','b','Linewidth',2)
    %     plot(x,PET7','Linewidth',2)
    %     legend('UB','LB','Mean','PETmean','Location','northwest')
    %     axis([70,150,0,1])
    %     box on
    %     grid on
    
    figure;
    hold on;
    title('Patient07')
    set(gca,'Fontsize',fs);
    jbfill(x,LB7',UB7',b4,a4,0,0.5);
    plot(x,PET7','Linewidth',2)
    plot(x,Mean7s','Linewidth',2)
    legend('2std','PET','Mean','Location','northwest')
    axis([70,150,0,1])
    box on
    grid on
    print('P07','-djpeg')
    print('P07','-dpdf')
    
end;

%% Plot one pet on top of each other
% PETreg1 = zeros(3,1);
% PETreg7 = zeros(3,1);
% i = 1;
% j = 1;
%
% PET1 = PET1./max(PET1(:));
% PET7 = PET7./max(PET7(:));
%
% for ix = 1:256
%     if(PET7(ix) > 0 )
%         PETreg7(i) = PET7(ix);
%         i = i + 1;
%     end;
%     if(PET1(ix) > 0 )
%         PETreg1(j) = PET1(ix);
%         j = j + 1;
%     end;
% end;
%
% l1 = length(PETreg1);
% l7 = length(PETreg7);
%
%
% figure;
% hold on
% set(gca,'Fontsize',fs);
% % plot(2:l1+1,PETreg1')
% % plot(1:l7,PETreg7')
%  plot(x,PET1')
%  plot(x,PET7')
%  axis([80,220,0,1])
% legend('01','07')
% box on
% grid on


