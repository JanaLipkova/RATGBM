%========================================
%
%  Read in output of CMA from outxrecentbests. file,
%  and rearange data so they can be plot with TMCMC ploting toos and with
%  Rscript, also rescale to correct parameter space
%
%========================================

%function CMAdataConversion

close all; clear all; clc

addpath('../../lib/jbfill')
addpath('../../lib/')

dataType = 'ALL';

if(dataType == 'PET')
    scaling = [1, 10, 100000, 1, 100, 100];
    fname = sprintf('P01_CMA_PET.txt');
    mydata = importdata('../../../CMA/Patients/Patient01/Patient01PET/outcmaesxrecentbest.dat');
else
    scaling = [1, 10, 100000, 1, 100, 100, 100, 100, 10];
    fname = sprintf('P01_CMA_All.txt');
    mydata = importdata('../../../CMA/Patients/Patient01/Patient01ALL/outcmaesxrecentbest.dat');
end;

[Nx,Ny] = size(mydata.data);

% 1) Rearange data : 5th column is likelihood, 6:end are parameters
outdata = zeros(Nx,Ny-4);
outdata(:,end) =  mydata.data(:,5).*(-1) ;

j=1;
for i=6:Ny
    outdata(:,j) = mydata.data(:,i);
    j=j+1;
end;

%2) remove zeros from the 1st row
outdata = outdata(2:end, :);

% 3) rescale to parameter space
for i = 1:Ny-5
    outdata(:,i) = exp( outdata(:,i) ) .* scaling(i);
end;

D = outdata(:,1);
ss = find(1==D(:))
stop = min(ss) - 1

 outdata = outdata(1:stop,:);

% % 3) save output for Rscript
fid = fopen(fname, 'wt'); % Open for writing
[Nx,Ny] = size(outdata);
for i=1:Nx
    for j=1:Ny
        fprintf(fid, '%d ', outdata(i,j));
    end
    fprintf(fid, '\n ');
end;
fclose(fid);

bestC = find( max(outdata(:,end)) == outdata(:,end));
best = outdata(bestC,:)
meanData = mean(outdata)
varData = var(outdata)
stdData = sqrt(varData)