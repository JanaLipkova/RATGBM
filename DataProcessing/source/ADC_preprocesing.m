% ===============================================
%  ADC Preprocessing
%
%   - ADC signal is inversly proportional to tumour cell density
%   - 1) subtruck healthy brain tissue signal (i.e. mean signal from healthy hemisphere)
%   - 2) Inverse
%   - 3) Normalised
%   - 4) restrict to tumour regions and vice-versa
%
%
%  Copyright: Jana Lipkova
%  jana.lipkova@tum.de
% ===============================================

% function ADC_preprocesing

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')



rID = 42;
day=9;
dataBase   = ['../../RAT_DATA_F98/M',num2str(rID)];

%0) Set path to data
dataPath = [dataBase,'/Modalities/day',num2str(day,'%02d')];
tumPath  = [dataBase,'/Tumor-ROI/corrected/'];
maskPath = [dataBase,'/Mask'];


%1) Read in tissue segmentaiotns
adc  = MRIread([dataPath,'/M',num2str(rID),'J',num2str(day,'%02d'),'-CoregADC.nii']);
t2w  = MRIread([tumPath, '/M',num2str(rID),'_J',num2str(day,'%02d'),'_T2w.nii.gz']);
mask = MRIread([maskPath,'/M',num2str(rID),'J',num2str(day,'%02d'),'-Mask.nii']);

% 2) Restrict
adc.vol = adc.vol.*mask.vol;
t2w_margin = t2w.vol;%removeNvoxelsFromBorder(t2w.vol,2,0);
adc.vol = adc.vol.* t2w_margin;

% 2) Inverse the ADC values
[Nx,Ny,Nz] = size(adc.vol);
iadc = zeros(Nx,Ny,Nz);

for iz=1:Nz
    for iy=1:Ny
        for ix=1:Nx
            if(adc.vol(ix,iy,iz)>0)
                iadc(ix,iy,iz) = 1./adc.vol(ix,iy,iz);
            end;
        end;
    end;
end;

% 3) normalised
niadc = iadc./max(iadc(:));
vi(niadc)

% 4) Apply gaussian smootheing filter
Iblur15 = imgaussfilt(niadc,1.5);
vi(Iblur15)

Iblur2 = imgaussfilt(niadc,2);
vi(Iblur2)

% 5) normallised again
Iblur15 = Iblur15./max(Iblur15(:));

% mask.vol = Iblur15;
% MRIwrite(mask,[dataPath,'/M',num2str(rID),'J',num2str(day,'%02d'),'-CoregInverseADC.nii']);
