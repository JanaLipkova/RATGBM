%=======================================
%
%  Compute P% confidence interval
%
%=======================================
% Input: 
%   x = 1D vector of data
%   P = percentage of confidence itnervals, i.e. 0.95
% Output:
%   xd = down value
%   xu = upper value
%   median, mean and std
%=======================================

function [xd,xu,x_median, x_mean, x_std] = confidenceIntervals(x, P)

xs=sort(x);
F=linspace(0,1,length(xs)); %Primitive function
clf;plot(xs,F,'.'); %Plot primitive distribution F=int_0^x p(x)dx
Fd=(1-P)/2; %Down value
Fu=1-Fd; %Upper value
%interpolate xu corr. to Pu etc.


xd=pchip(F,xs,Fd);
xu=pchip(F,xs,Fu);
x_median=pchip(F,xs,0.5);
xm=mean(x);
dx=std(x);

bounds = [xd,xu]

figure;
hold on
title([num2str(P*100),'% CI=[',num2str(xd),' , ',num2str(xu),'], median=',num2str(x_median)])
hist(xs,30);hold on;
plot([xd xd],[0 250],'r',...
[xu xu],[0 50],'r',[xm xm],[0 100],'g') 