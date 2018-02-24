function KDEtesting3

addpath('/Users/lipkova/WORK/MATLAB_codes/kde2d');
addpath('/Users/lipkova/WORK/MATLAB_codes');

% kde2 solver
% generate a Gaussian mixture with distant modes
data=[randn(100,1), randn(100,1)/4;
    randn(100,1)+18, randn(100,1);
    randn(100,1)+15, randn(100,1)/2-18;];
% call the routine
[bandwidth,density,X,Y]=kde2d(data);

% plot the data and the density estimate
fid=1
figure(fid)
surf(X,Y,density,'LineStyle','none'), view([0,60])
colormap hot, hold on, alpha(.8)
set(gca, 'color', 'blue');
plot(data(:,1),data(:,2),'w.','MarkerSize',5)
fid=fid+1;

size(data)
% original kde
p = kde( data', [0.3;0.3]);
figure(fid)
plot(p);

f=figure(fid)
pcolor(hist(p));
shading flat
fid = fid+1;
saveas(f,'tsts.pdf')

figure(fid)
mesh(hist(p));