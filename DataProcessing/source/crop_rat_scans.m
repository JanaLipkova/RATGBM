%===============================================
%
%   Crop data to save computation costs
%
%
%  Copyright: Jana Lipkova
%  jana.lipkova@tum.de
%===============================================


% function crop_rat_scans

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')

rID = 30 %42 38 36 34 30 29
day_list = [9,11,14,16];
bAnatomy = 1;
bRATS=1;
bATLAS=0;


if(bRATS)
    for day = day_list
        day
        % 0) set up path to data
        inputRatPath  = ['../../RAT_DATA_F98/M',num2str(rID),'/isotropic/'];
        inputModPath  = [inputRatPath,'modalities/'];
        inputAnatPath = [inputRatPath,'anatomy/'];
        inputTumPath  = [inputRatPath,'tumour/'];
        outputPath    = ['../../RAT_DATA_F98/M',num2str(rID),'/isotropic_cropped/'];
        
        % 1) read in data
        T2wm = MRIread([inputModPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-Coreg_Anat-masked-iso.nii'] );
        DCE  = MRIread([inputModPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-CoregDCE-AUC-iso.nii'] );
        
        mask  = MRIread([inputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-Mask-iso.nii'] );
        T2ws  = MRIread([inputTumPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-T2w-iso.nii'] );
        DCEs  = MRIread([inputTumPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-DCE-iso.nii'] );
        
        DCE.vol = DCE.vol.*mask.vol;
        
        
        % 2) crop
        mask_tmp = mask.vol;
        T2wm.vol = cropData(mask_tmp, T2wm.vol);
        DCE.vol  = cropData(mask_tmp, DCE.vol);
        mask.vol = cropData(mask_tmp, mask.vol);
        T2ws.vol = cropData(mask_tmp, T2ws.vol);
        DCEs.vol = cropData(mask_tmp, DCEs.vol);
        
        
        % 3) Add 4 ghosts points around the cropped volume
        T2wm.vol = padarray(T2wm.vol,[4 4 4],0,'both');
        DCE.vol  = padarray(DCE.vol, [4 4 4],0,'both');
        mask.vol = padarray(mask.vol,[4 4 4],0,'both');
        T2ws.vol = padarray(T2ws.vol,[4 4 4],0,'both');
        DCEs.vol = padarray(DCEs.vol,[4 4 4],0,'both');
        
        size(T2wm.vol)
        
        % 5) save cropped data to nii
        MRIwrite(T2wm, [outputPath,'nifty/modalities/M',num2str(rID),'J',num2str(day,'%02d'),'-Coreg_Anat-masked-iso-crop.nii']);
        MRIwrite(DCE,  [outputPath,'nifty/modalities/M',num2str(rID),'J',num2str(day,'%02d'),'-CoregDCE-AUC-iso-crop.nii']);
        MRIwrite(mask,  [outputPath,'nifty/anatomy/M',num2str(rID),'J',num2str(day,'%02d'),'-Mask-iso-crop.nii']);
        MRIwrite(T2ws,  [outputPath,'nifty/tumour/M',num2str(rID),'J',num2str(day,'%02d'),'-T2w-iso-crop.nii']);
        MRIwrite(DCEs,  [outputPath,'nifty/tumour/M',num2str(rID),'J',num2str(day,'%02d'),'-DCE-iso-crop.nii']);
        
        
        % 6) save cropped data into binary
        typeId = 1; % 0 double, 1 float
        writeMatrix(T2wm.vol,  [outputPath,'binary/modalities/M',num2str(rID),'J',num2str(day,'%02d'),'-Coreg_Anat-masked-iso-crop.dat'],typeId)
        writeMatrix(DCE.vol,   [outputPath,'binary/modalities/M',num2str(rID),'J',num2str(day,'%02d'),'-CoregDCE-AUC-iso-crop.dat'],typeId);
        writeMatrix(mask.vol,  [outputPath,'binary/anatomy/M',num2str(rID),'J',num2str(day,'%02d'),'-Mask-iso-crop.dat'],typeId);
        writeMatrix(T2ws.vol,  [outputPath,'binary/tumour/M',num2str(rID),'J',num2str(day,'%02d'),'-T2w-iso-crop.dat'],typeId);
        writeMatrix(DCEs.vol,  [outputPath,'binary/tumour/M',num2str(rID),'J',num2str(day,'%02d'),'-DCE-iso-crop.dat'],typeId);
        
    end;
    
    
    % Apply the same for the anatomy
    if(bAnatomy)
        day=9;
        csf = MRIread([inputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-csf-iso.nii'] );
        wm  = MRIread([inputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-wm-iso.nii'] );
        gm  = MRIread([inputAnatPath,'M',num2str(rID),'J',num2str(day,'%02d'),'-gm-iso.nii'] );
        
        csf.vol = cropData(mask_tmp, csf.vol);
        wm.vol  = cropData(mask_tmp, wm.vol);
        gm.vol  = cropData(mask_tmp, gm.vol);
        
        csf.vol = padarray(csf.vol,[4 4 4],0,'both');
        wm.vol  = padarray(wm.vol,[4 4 4],0,'both');
        gm.vol  = padarray(gm.vol,[4 4 4],0,'both');
        
        MRIwrite(csf, [outputPath,'nifty/anatomy/M',num2str(rID),'J',num2str(day,'%02d'),'-csf-iso-crop.nii']);
        MRIwrite(wm,  [outputPath,'nifty/anatomy/M',num2str(rID),'J',num2str(day,'%02d'),'-wm-iso-crop.nii']);
        MRIwrite(gm,  [outputPath,'nifty/anatomy/M',num2str(rID),'J',num2str(day,'%02d'),'-gm-iso-crop.nii']);
        
        typeId=1;
        writeMatrix(csf.vol,  [outputPath,'binary/anatomy/M',num2str(rID),'J',num2str(day,'%02d'),'-csf-iso-crop.dat'],typeId)
        writeMatrix(wm.vol,  [outputPath,'binary/anatomy/M',num2str(rID),'J',num2str(day,'%02d'),'-wm-iso-crop.dat'],typeId)
        writeMatrix(gm.vol,  [outputPath,'binary/anatomy/M',num2str(rID),'J',num2str(day,'%02d'),'-gm-iso-crop.dat'],typeId)
    end;
    
end;



if(bATLAS)
    inputPath = '../../AtlasRat/atlas/';
    atlas = MRIread([inputPath,'rat_atlas.nii']);
    mask  = MRIread([inputPath,'rat_atlas_mask.nii']);
    wm    = MRIread([inputPath,'rat_atlas_wm.nii']);
    gm    = MRIread([inputPath,'rat_atlas_gm.nii']);
    csf   = MRIread([inputPath,'rat_atlas_csf.nii']);
    
    %cut in the z-direction to save performance
    mask_tmp = mask.vol;
    mask_tmp(:,:,1:40)=0;
    
    atlas.vol = cropData(mask_tmp, atlas.vol);
    mask.vol  = cropData(mask_tmp, mask.vol);
    wm.vol    = cropData(mask_tmp, wm.vol);
    gm.vol    = cropData(mask_tmp, gm.vol);
    csf.vol   = cropData(mask_tmp, csf.vol);
    
    % pad with 4 zeros ghost points
    atlas.vol = padarray(atlas.vol,[4 4 4],0,'both');
    mask.vol = padarray(mask.vol,[4 4 4],0,'both');
    wm.vol = padarray(wm.vol,[4 4 4],0,'both');
    gm.vol = padarray(gm.vol,[4 4 4],0,'both');
    csf.vol = padarray(csf.vol,[4 4 4],0,'both');
    
    % make segmentations binary
    LB=0.3;
    mask.vol(mask.vol(:)<=LB) = 0;
    mask.vol(mask.vol(:)> LB) = 1;
    wm.vol(wm.vol(:)<=LB) = 0;
    wm.vol(wm.vol(:)> LB) = 1;
    gm.vol(gm.vol(:)<=LB) = 0;
    gm.vol(gm.vol(:)> LB) = 1;
    csf.vol(csf.vol(:)<=LB) = 0;
    csf.vol(csf.vol(:)> LB) = 1;
    
    MRIwrite(atlas,[inputPath,'rat_atlas_crop.nii']);
    MRIwrite(mask, [inputPath,'rat_atlas_mask_crop.nii']);
    MRIwrite(wm,   [inputPath,'rat_atlas_wm_crop.nii']);
    MRIwrite(gm,   [inputPath,'rat_atlas_gm_crop.nii']);
    MRIwrite(csf,  [inputPath,'rat_atlas_csf_crop.nii']);
    
    typeId=1; % 0-double, 1 float
    writeMatrix(atlas.vol, [inputPath,'rat_atlas_crop.dat'],typeId)
    writeMatrix(mask.vol,  [inputPath,'rat_atlas_mask_crop.dat'],typeId)
    writeMatrix(wm.vol,    [inputPath,'rat_atlas_wm_crop.dat'],typeId)
    writeMatrix(gm.vol,    [inputPath,'rat_atlas_gm_crop.dat'],typeId)
    writeMatrix(csf.vol,   [inputPath,'rat_atlas_csf_crop.dat'],typeId)
    
end;
