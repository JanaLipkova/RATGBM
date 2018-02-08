%===================================
%
%   Test Synthetic data
%
%==================================



% function testDataCheck

% close all;
clc


addpath('/Users/lipkova/WORK/UQ/UQ_VillaGarbald/Matlab2C/matrixMatlab2Cpp/matlab/')


% A=loadMatrix('/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient01/FET1_1b.dat');
A=loadMatrix('SyntheticData/HGG_data.dat');
size(A)

m=0;
c=1;
siceID=64
v=[1 2 3 4 5 6];
for step = v
    
    name = ['Sphere_',num2str(1),'.dat']
    S=loadMatrix(name);
    
    B=zeros(size(A));
    
    S = S';
    out = zeros(6,2);
    Npoints = 0;
    
    for j = 1:size(S,1)
        ix = S(j,1);
        iy = S(j,2);
        iz = S(j,3);
        
        B(ix,iy,iz) = 10;
        
        if( (mod(ix,step)==0) && (mod(iy,step)==0) && (mod(iz,step)==0) )
        Npoints = Npoints + 1;
        end;
        
        if((iz==siceID)&&(mod(ix,step)==0)&&(mod(iy,step)==0))
            m=m+1;
            out(m,1) = ix;
            out(m,2) = iy;
        end;
        
    end;
 
    
    figure,
    hold on
    pcolor(A(:,:,siceID))
    scatter(out(:,2), out(:,1),'*w')
    name2=['Step=',num2str(step),' Points=',num2str(Npoints),' Slice=',num2str(siceID)];
    title(name2)
    colorbar;
    colormap jet
    saveas(gca,[name2,'.jpg'])

end;


%%
figure,
pcolor(A(:,:,64))
colorbar
colormap jet
saveas(gca,'tumor.jpg')


% figure(fid)
% pcolor(A(:,:,63))
% fid = fid+1;
%
% figure(fid)
% pcolor(A(:,:,62))
% fid = fid+1;
%
% figure(fid)
% pcolor(A(:,:,61))
% fid = fid+1;



