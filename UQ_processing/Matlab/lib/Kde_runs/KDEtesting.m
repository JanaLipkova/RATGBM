function KDEtesting

%addpath('/Users/lipkova/Documents/MATLAB/jbfill')
addpath('/Users/lipkova\ 1/Documents/MATLAB/jbfill')
addpath('../')

% kde(points,ksize[,weights][,type]) -- nonpara'metric density estimate of a pdf
%
%  points is the [Ndim x Npoints] array of kernel locations
%  ksize may be a scalar, [Ndim x 1], [Ndim x Npoints], or
%    a string (for data-based methods; see @kde/ksize for allowed methods)
%  weights is [1 x Npoints] and need not be pre-normalized
%  type can be one of: 'Gaussian', 'Laplacian', 'Epanetchnikov'
%    (only 1st letter required) (Gaussian by default)
%


close all;
fid=1;
N=1000;
x = randn(2,N)

p = kde( x, [0.6;0.6]);

figure(fid)
plot(p);
fid = fid+1;

figure(fid)
mesh(hist(p));
fid = fid+1;

% p = ksize(p, 'hall'); % Plug-in type estimator (estimates each dim. separately)

figure(fid)
mesh(hist(p));
fid = fid+1;

figure(fid)
pcolor(hist(p));
shading flat
fid = fid+1;

xl = min(x(1,:))
xu = max(x(1,:))
X=linspace(xl,xu,200);

yl = min(x(2,:))
yu = max(x(2,:))
Y=linspace(yl,yu,200);

figure(fid)
pcolor(X,Y,hist(p));
shading flat
fid = fid+1;

size(hist(p))
% statistic
m1 = mean(p)
C=covar(p)

%% marginalisation
p2 = kde(x(1,:), 0.3);
size(p2)

figure(fid)
g=plot(p2);
set(g,'Color',[0.6,0.0,1.0],'Linewidth',2)
fid = fid+1;
double(p2)

pts = getPoints(p2);
N = 200;   range = [min(pts),max(pts)];
bins=linspace(range(1),range(2),N);
y = evaluate(p2,bins);

figure(fid)
plot(bins,y);
fid=fid+1;

z=zeros(size(y));
figure(fid)
[fillhandle,msg]=jbfill(bins,y,z,rand(1,3),rand(1,3),0,rand(1,1))




