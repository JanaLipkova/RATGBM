% function OverlaySamples



addpath('../../lib/jbfill')
addpath('../../lib/tightfig')
addpath('../../lib/')
addpath('../../lib/altmany-export_fig-4c015d5')

%% 1) Path to cirgen_db file + index of generation to be plot
%-------------------------------------------------------------
pathname = '/../../../../Volumes/brutus_scratch/TMCMC/PatientCases/';

%-----------------------
Patient = 1

GenId = 14;
myfilename = sprintf('Patient01/V2/P01_V2_4K/curgen_db_%03d.txt',GenId);

myfilename = [pathname,myfilename]
mydata = importdata(myfilename);

mydata(:,1) = exp( mydata(:,1) ) * 100;  % convert from cm^2 to mm^2
mydata(:,2) = exp( mydata(:,2) );

dataP01 = mydata(:,1:2);
clear mydata;
%-----------------------

Patient = 7
GenId = 16;
myfilename = sprintf('Patient07/V2/P07_V2_4K/curgen_db_%03d.txt',GenId);

myfilename = [pathname,myfilename];
mydata = importdata(myfilename);

mydata(:,1) = exp( mydata(:,1) ) * 100;  % convert from cm^2 to mm^2
mydata(:,2) = exp( mydata(:,2) );

dataP07 = mydata(:,1:2);
clear mydata;
%-----------------------

Patient = 11
GenId = 21;
myfilename = sprintf('Patient11/V2/P11_V2_4K_b/curgen_db_%03d.txt',GenId);

myfilename = [pathname,myfilename];
mydata = importdata(myfilename);

mydata(:,1) = exp( mydata(:,1) ) * 100;  % convert from cm^2 to mm^2
mydata(:,2) = exp( mydata(:,2) );

dataP11 = mydata(:,1:2);
clear mydata;
%-----------------------

Patient = 22
GenId = 19;
myfilename = sprintf('Patient22/V2/P22_V2_4K_b/curgen_db_%03d.txt',GenId);

myfilename = [pathname,myfilename];
mydata = importdata(myfilename);

mydata(:,1) = exp( mydata(:,1) ) * 100;  % convert from cm^2 to mm^2
mydata(:,2) = exp( mydata(:,2) );

dataP22 = mydata(:,1:2);
clear mydata;
%-----------------------



% Synthetic = 1
% GenId = 22
% myfilename = sprintf('../../../../../Volumes/brutus_scratch/TMCMC/SpaceSearching/SynthethicAll_All/SynthethicAll_All_4K_Lustre_C/curgen_db_%03d.txt',GenId);
% mydata = importdata(myfilename);
% 
% mydata(:,1) = exp( mydata(:,1) ) * 100;  % convert from cm^2 to mm^2
% mydata(:,2) = exp( mydata(:,2) );
% 
% dataS1 = mydata(:,1:2);
% clear mydata;
% %-----------------------
% 
% Synthetic = 2
% GenId = 28;
% myfilename = sprintf('SyntheticBig/SynthethicAll_Big_4K_Lustre_NewB/curgen_db_%03d.txt',GenId);
% 
% myfilename = [pathname,myfilename];
% mydata = importdata(myfilename);
% 
% mydata(:,1) = exp( mydata(:,1) ) * 100;  % convert from cm^2 to mm^2
% mydata(:,2) = exp( mydata(:,2) );
% 
% dataS2 = mydata(:,1:2);
% clear mydata;
% %-----------------------



%% Plot HGG space

% DR = [1,1000] / 365;
rhoR = [0.1,100] / 365;
DR = [ 1/365    3.800];



% Points for lines in plot
p1x = [2,1000];
p1y = [0.1,50];

p2x = [1,200];
p2y = [0.5,100];
p3x = [10, 1];
p3y = [0.1,1];
p4x = [250, 1];
p4y = [0.1,25];
p5x = [1000, 100];
p5y = [10,100];
p6x = [5  ,1 ];
p6y = [0.1,0.5];
p7x = [20  ,1 ];
p7y = [0.1, 2];

% plotting set up
close all

fs=30;
lw = 2;
pointsize=25;



%color maps
%purple
c1 = [238,130,238] / 255 ; 
c2 = [186, 85, 211] /255; 
c3 = [160, 32, 240] /255; 

% c4 = [1.0, 0.2, 0.0] ; % orange
% c5 = [1.0, 0.6, 0.0] ; % blue
% c6 = [1.0, 0.6, 0.2] ; % blue

c4 = [0.0, 0.6, 0.2] ; % purple
c5 = [0.0, 1.0, 0.2] ; % blue
c6 = [0.2, 0.4, 0.0] ; % blue

%%
figure();
set(gca,'Fontsize',fs);
hold on
plot(p1x / 365,p1y / 365,'k-', 'Linewidth',lw)
plot(p2x / 365,p2y / 365,'k-','Linewidth',lw)
% plot(p3x / 365,p3y / 365,'k--','Linewidth',lw)
plot(p4x / 365,p4y / 365,'k--','Linewidth',lw)
plot(p5x / 365,p5y / 365,'k--','Linewidth',lw)
plot(p6x / 365,p6y / 365,'k--','Linewidth',lw)
plot(p7x / 365,p7y / 365,'k--','Linewidth',lw)

scatter(dataP01(:,1), dataP01(:,2),pointsize,c1,'filled' );
scatter(dataP07(:,1), dataP07(:,2),pointsize,c2,'filled');
% scatter( dataS1(:,1),  dataS1(:,2),pointsize,c3,'filled');


scatter(dataP11(:,1), dataP11(:,2),pointsize,'b','filled');
scatter(dataP22(:,1), dataP22(:,2),pointsize,c5,'filled');
% scatter( dataS2(:,1),  dataS2(:,2),pointsize,c6,'filled');



axis([DR(1) ,DR(2),rhoR(1), rhoR(2)])
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'linewidth',3)
axis square;

xlabel('D(mm^2\day)');
ylabel('\rho (1\day)')
box on;  grid on


set(gca, 'Color', 'None');
set(gcf, 'Position', [100, 100, 1049, 895]);
set(gcf,'papersize',[400,400]);

export_fig testJana.png -transparent;







