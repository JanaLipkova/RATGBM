%===============================================
%
%   Make data isotropic
%
%   - resize data into uniform resolution
%
%  Copyright: Jana Lipkova
%  jana.lipkova@tum.de
%===============================================


% function make_data_isotropic

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')

% 0) set input paht
rID = 38; %[30,34,36,38,42]   %42,38,36,34,30,29
bAnatomy = 1;
day_list =  [9,11,14,16];

for day = day_list
    day
    
    inputRatPath   = ['../../RAT_DATA_F98/M',num2str(rID),'/'];
    inputModPath   = [inputRatPath,'Modalities/day',num2str(day,'%02d'),'/'];
    inputAnatPath  = [inputRatPath,'Anat-VOI/'];
    inputMaskPath  = [inputRatPath,'Mask/'];
    inputTumPath   = [inputRatPath,'Tumor-ROI/corrected/new2/'];
    
    outputTumPath  = [inputRatPath,'isotropic/tumour/'];
    outputAnatPath = [inputRatPath,'isotropic/anatomy/'];
    outputModPath  = [inputRatPath,'isotropic/modalities/'];
    
    % 1 read in modalities, mask, tumour    
    T2wm = MRIread([inputModPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-Coreg01_Anat-masked.nii']);
%     T2w  = MRIread([inputModPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-Coreg01_Anat.nii']);
    DCE  = MRIread([inputModPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-CoregDCE-AUC.nii']);
    mask = MRIread([inputMaskPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-Mask.nii']);
    T2ws = MRIread([inputTumPath,'M',num2str(rID),'_J',num2str(day,'%02d'),'_T2w.nii.gz']);
    DCEs = MRIread([inputTumPath,'M',num2str(rID),'_J',num2str(day,'%02d'),'_DCE.nii.gz']);
    
    % 2) make data isotropic
    [Nx,Ny,Nz] = size(T2wm.vol);
    Nz = Nz*6;
    newRes = [Nx,Ny,Nz];
    interp = '*linear';   %'linear', 'nearest', 'cubic', 'makima', or 'spline'
    
    T2wm.vol = resize3d(T2wm.vol, newRes,interp);
%     T2w.vol  = resize3d(T2w.vol,  newRes,interp);
    DCE.vol  = resize3d(DCE.vol,  newRes,interp);
    mask.vol = resize3d(mask.vol, newRes,interp);
    T2ws.vol = resize3d(T2ws.vol, newRes,interp);
    DCEs.vol = resize3d(DCEs.vol, newRes,interp);
    
    % 3) Threashold binary masks, restrict inside the mask
    LB=0.5;
    mask.vol(mask.vol(:)<=LB) = 0;
    mask.vol(mask.vol(:)> LB) = 1;
    
    T2ws.vol = T2ws.vol.*mask.vol;
    DCEs.vol = DCEs.vol.*mask.vol;
    
    T2ws.vol(T2ws.vol(:)<=LB) = 0;
    T2ws.vol(T2ws.vol(:)> LB) = 1;
    
    DCEs.vol(DCEs.vol(:)<=LB) = 0;
    DCEs.vol(DCEs.vol(:)> LB) = 1;
    
    % 4) save
    MRIwrite(T2wm, [outputModPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-Coreg_Anat-masked-iso.nii.gz']);    
%     MRIwrite(T2w,  [outputModPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-Coreg_Anat-iso.nii.gz']);
    MRIwrite(DCE,  [outputModPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-CoregDCE-AUC-iso.nii.gz']);
    
    MRIwrite(mask, [outputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-Mask-iso.nii.gz']);
    MRIwrite(T2ws, [outputTumPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-T2w-iso.nii.gz']);
    MRIwrite(DCEs, [outputTumPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-DCE-iso.nii.gz']);
end;

%% Anatomies
if (bAnatomy)
    day=9;
    
    csf = MRIread([inputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-VOI-Atlas_CSF.nii']);
    wm  = MRIread([inputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-VOI-Atlas_WM.nii']);
    gm  = MRIread([inputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-VOI-Atlas_GM.nii']);
    
    csf.vol = resize3d(csf.vol, newRes,interp);
    wm.vol  = resize3d(wm.vol, newRes,interp);
    gm.vol  = resize3d(gm.vol, newRes,interp);
    
    
    csf.vol = csf.vol .* mask.vol;
    wm.vol  = wm.vol  .* mask.vol;
    gm.vol  = gm.vol  .* mask.vol;
    
    csf.vol(csf.vol(:)<=LB) = 0;
    csf.vol(csf.vol(:)> LB) = 1;
    
    wm.vol(wm.vol(:)<=LB) = 0;
    wm.vol(wm.vol(:)> LB) = 1;
    
    gm.vol(gm.vol(:)<=LB) = 0;
    gm.vol(gm.vol(:)> LB) = 1;
    
    % correct anatomy
    [Nx,Ny,Nz] = size(gm.vol);
    for iz=1:Nz
        for iy=1:Ny
            for ix=1:Nx
                if(csf.vol(ix,iy,iz) > 0.5)
                    gm.vol(ix,iy,iz) = 0;
                    wm.vol(ix,iy,iz) = 0;
                else
                    if(mask.vol(ix,iy,iz) > 0.5)
                        p_wm = wm.vol(ix,iy,iz) / (wm.vol(ix,iy,iz) + gm.vol(ix,iy,iz));
                        p_gm = gm.vol(ix,iy,iz) / (wm.vol(ix,iy,iz) + gm.vol(ix,iy,iz));
                        wm.vol(ix,iy,iz) = p_wm;
                        gm.vol(ix,iy,iz) = p_gm;
                    end;
                end;
            end;
        end;
    end;
    
    MRIwrite(csf, [outputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-csf-iso.nii.gz']);
    MRIwrite(wm,  [outputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-wm-iso.nii.gz']);
    MRIwrite(gm,  [outputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-gm-iso.nii.gz']);
end;




