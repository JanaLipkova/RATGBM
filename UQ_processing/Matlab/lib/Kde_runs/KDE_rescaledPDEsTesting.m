%function KDE_rescaledPDEsTesting

addpath('/Users/lipkova/WORK/MATLAB_codes');
addpath('/Users/lipkova/Documents/MATLAB/jbfill')

randn('state',1);
close all;

% KDE set ups
N = 200;
M = 1000;

% generate a Gaussian mixture with distant modes

data=[randn(M,1), randn(M,1)];

fid=1;
toBeSafedAs ='Original';
plotKDE(data,fid,N,toBeSafedAs)

fid=2;
% data=data*10;
% toBeSafedAs = 'Rescled10'
data=data.^1/3;
toBeSafedAs = 'RescledSqrt3'
plotKDE(data,fid,N,toBeSafedAs)




