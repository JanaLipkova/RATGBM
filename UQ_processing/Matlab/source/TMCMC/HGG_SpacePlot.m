% function HGG_SpacePlot
close all

%% [mm/year], [1/year]
DR = [1,1000];
rhoR = [0.1,100];

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


%plotting set up
fs=30;
lw = 2;


figure;
set(gca,'Fontsize',fs);
hold on
plot(p1x,p1y,'k-', 'Linewidth',lw)
plot(p2x,p2y,'k-','Linewidth',lw)
% plot(p3x,p3y,'k--','Linewidth',lw)
plot(p4x,p4y,'k--','Linewidth',lw)
plot(p5x,p5y,'k--','Linewidth',lw)
plot(p6x,p6y,'k--','Linewidth',lw)
plot(p7x,p7y,'k--','Linewidth',lw)


axis([DR(1),DR(2),rhoR(1), rhoR(2)])
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'linewidth',3) 
axis square;

xlabel('D[mm^2\yr]');
ylabel('rho [1\yr]');
box on;  grid on

set(gcf, 'Position', [100, 100, 1049, 895]);
set(gcf,'papersize',[300,300]);
set(gcf,'PaperPositionMode','auto')

% print('HGG_mmYear2','-djpeg')



%% [mm/day] [1/day]
% 
DR = [1,1000] / 365;
rhoR = [0.1,100] / 365;

figure;
set(gca,'Fontsize',fs);
hold on
plot(p1x / 365,p1y / 365,'k-', 'Linewidth',lw)
plot(p2x / 365,p2y / 365,'k-','Linewidth',lw)
% plot(p3x / 365,p3y / 365,'k--','Linewidth',lw)
plot(p4x / 365,p4y / 365,'k--','Linewidth',lw)
plot(p5x / 365,p5y / 365,'k--','Linewidth',lw)
plot(p6x / 365,p6y / 365,'k--','Linewidth',lw)
plot(p7x / 365,p7y / 365,'k--','Linewidth',lw)


axis([DR(1) ,DR(2),rhoR(1), rhoR(2)])
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'linewidth',3) 
axis square;

xlabel('D(mm^2\day)');
ylabel('\rho (1\day)')
box on;  grid on

set(gcf, 'Position', [100, 100, 1049, 895]);
set(gcf,'papersize',[400,400]);
set(gcf,'PaperPositionMode','auto')
% print('HGG_mmDay2','-djpeg')
