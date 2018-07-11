% ===============================================
%  Create brain mask
%
%   - used masked T2w scans to create brain masks
%   AND
%   - original tissue segmentations to create brain mask
%
%  Copyright: Jana Lipkova
%  jana.lipkova@tum.de
% ===============================================

% function create_brain_masks

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')


%0) select case and set-up paths
rList = [29,30,34,36,38,42];

for rID = rList
    
    rID
    dataBase   = '../../RAT_DATA_F98/';
    bModality = 0; % 0 for anatomy, 1 for modality
    
    if(bModality)
        
        modalitiesPath = [dataBase,'M',num2str(rID),'/Modalities/'];
        outMasksPath   = [dataBase,'M',num2str(rID),'/Mask/'];
        
        % 0) Read in masked T2w image from each time point
        mod_09 = MRIread([modalitiesPath,'day09/M',num2str(rID),'J09-Coreg01_Anat-masked.nii']);
        mod_11 = MRIread([modalitiesPath,'day11/M',num2str(rID),'J11-Coreg01_Anat-masked.nii']);
        mod_14 = MRIread([modalitiesPath,'day14/M',num2str(rID),'J14-Coreg01_Anat-masked.nii']);
        mod_16 = MRIread([modalitiesPath,'day16/M',num2str(rID),'J16-Coreg01_Anat-masked.nii']);
        
        % 1) create brain masks
        mask_09 = mod_09.vol;
        mask_11 = mod_11.vol;
        mask_14 = mod_14.vol;
        mask_16 = mod_16.vol;
        
        cut = 0.9;
        
        mask_09(mask_09 > cut) = 1;
        mask_11(mask_11 > cut) = 1;
        mask_14(mask_14 > cut) = 1;
        mask_16(mask_16 > cut) = 1;
        
        mask_09(mask_09 < cut) = 0;
        mask_11(mask_11 < cut) = 0;
        mask_14(mask_14 < cut) = 0;
        mask_16(mask_16 < cut) = 0;
        
        
        % 2) save the masks
        mod_09.vol = mask_09;
        mod_11.vol = mask_11;
        mod_14.vol = mask_14;
        mod_16.vol = mask_16;
        
        MRIwrite(mod_09,[outMasksPath,'M',num2str(rID),'J09-Mask.nii']);
        MRIwrite(mod_11,[outMasksPath,'M',num2str(rID),'J11-Mask.nii']);
        MRIwrite(mod_14,[outMasksPath,'M',num2str(rID),'J14-Mask.nii']);
        MRIwrite(mod_16,[outMasksPath,'M',num2str(rID),'J16-Mask.nii']);
        
    else
        
        %0) Set path to data
        anatomyPath = [dataBase,'M',num2str(rID),'/Anat-VOI/'];
        
        %1) Read in tissue segmentaiotns
        gm  = MRIread([anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_GM.nii']);
        wm  = MRIread([anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_WM.nii']);
        csf = MRIread([anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_CSF.nii']);
        
        
        %2) Combined tissue into one brain mask
        mask = gm.vol + wm.vol + csf.vol;
        
        %3) Make the mask binary + save
        mask(mask(:)>0.5) = 1;
        
        %4) Save ouptput
        wm.vol = mask;
        MRIwrite(wm, [anatomyPath,'M',num2str(rID),'J09-Atlas_Mask.nii']);
        
    end;
end;