function plotKDE(data,fid,N,toBeSafedAs)


%plotting set ups
font = 12;
alpha = 0.5;
v1 = [0.4, 0.8, 1.0]; % light blue
v12 = [0.2, 0.2, 1.0]; % dark blue
v2 = [0.6, 0.6, 1.0]; % light purple
v22 = [0.6, 0.0, 0.8]; % dark purple


%grid values
xl = min(data(:,1));
xu = max(data(:,1));
yl = min(data(:,2));
yu = max(data(:,2));
X=linspace(xl,xu,N);
Y=linspace(yl,yu,N);
pointsize = 2;
bw1 = (xu - xl)/20
bw2 = (yu - yl)/20

% original kde
p = kde( data', [bw1;bw2]);
mean_value = mean(p)
Covariance = covar(p)

% plotting resuts
m = ksize(p, 'rot'); % "Rule of Thumb"; Silverman '86 / Scott '92
h=(hist(m,N,[2,1],[xl,xu;yl yu]));

f=figure(fid),
hold on;

subplot(4,2,1)
set(gca,'Fontsize',font);
hold on
scatter(data(:,1), data(:,2),pointsize,'*','Linewidth',3);
title('Samples')
axis([xl,xu,yl,yu]);
grid on;


subplot(4,2,2)
set(gca,'Fontsize',font);
pcolor(X,Y,h);
title('KDE estimate')
shading flat
colorbar

subplot(4,2,3)
set(gca,'Fontsize',font);
hold on
p1 = kde(data(:,1)', bw1);
% p1 = kde(h(1,:), bw1);
range = [xl,xu];
bins=linspace(range(1),range(2),N);
y = evaluate(p1,bins);
z=zeros(size(y));
[fillhandle,msg]=jbfill(bins,y,z,v1,v12,0,alpha);
set(fillhandle,'Linewidth',2)
title('Marginal p1')
xlim([xl xu])
grid on;


subplot(4,2,4)
set(gca,'Fontsize',font);
hold on
p2 = kde(data(:,2)', bw2);
% p2 = kde(h(2,:), bw2);
range = [yl,yu];
bins=linspace(range(1),range(2),N);
y = evaluate(p2,bins);
z=zeros(size(y));
[fillhandle,msg]=jbfill(bins,y,z,v2,v22,0,alpha);
set(fillhandle,'Linewidth',2)
title('Marginal p2')
xlim([yl yu])
grid on;


subplot(4,2,5)
set(gca,'Fontsize',font);
hold on
p1 = kde(h(1,:), bw1);
range = [xl,xu];
bins=linspace(range(1),range(2),N);
y = evaluate(p1,bins);
z=zeros(size(y));
[fillhandle,msg]=jbfill(bins,y,z,v1,v12,0,alpha);
set(fillhandle,'Linewidth',2)
title('Marginal KDE appr. p1')
xlim([xl xu])
grid on;


subplot(4,2,6)
set(gca,'Fontsize',font);
hold on
p2 = kde(h(2,:), bw2);
range = [yl,yu];
bins=linspace(range(1),range(2),N);
y = evaluate(p2,bins);
z=zeros(size(y));
[fillhandle,msg]=jbfill(bins,y,z,v2,v22,0,alpha);
set(fillhandle,'Linewidth',2)
title('Marginal KDE appro p2')
xlim([yl yu])
grid on;

t = [toBeSafedAs,'.jpg'];
saveas(f,t)