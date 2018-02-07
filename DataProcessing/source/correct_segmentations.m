%===============================================
%
%   Correct segmentations
%
%   - read in tumour segmentations (DCE and T2w) tumour
%   - restrict segmentations only inside the brain mask
%   - exclude segmentations from CSF regions
%   - apply time inclusion, i.e. if voxel is segmented in time 1, should
%   be segmented as a tumour in the next time point
%   - apply modalities inclusion, if segmented as lesion in DCE, should be
%   segmented inside T2w as well
%
%
%  Copyright: Jana Lipkova
%  jana.lipkova@tum.de
%===============================================


%function correct_segmentations

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')


%0) select case and set-up paths
rID = 34;
dataBase    = ['../../RAT_DATA_F98/M',num2str(rID),'/'];
MasksPath   = [dataBase, 'Mask/'];
tissuePath  = [dataBase, 'Anat-VOI/'];
tumourPath  = [dataBase, 'Tumor-ROI/'];


% 1) Read in masks, csf, tumour segm.
mask_09 = MRIread([MasksPath,'M',num2str(rID),'J09-Mask.nii']);
mask_11 = MRIread([MasksPath,'M',num2str(rID),'J11-Mask.nii']);
mask_14 = MRIread([MasksPath,'M',num2str(rID),'J14-Mask.nii']);
mask_16 = MRIread([MasksPath,'M',num2str(rID),'J16-Mask.nii']);

% csf     = MRIread([tissuePath,'M',num2str(rID),'J08-VOI-Atlas_CSF.nii']);

T2w_09  = MRIread([tumourPath,'M',num2str(rID),'_J09_Tumeur-ROI.nii']);
T2w_11  = MRIread([tumourPath,'M',num2str(rID),'_J11_Tumeur-ROI.nii']);
T2w_14  = MRIread([tumourPath,'M',num2str(rID),'_J14_Tumeur-ROI.nii']);
T2w_16  = MRIread([tumourPath,'M',num2str(rID),'_J16_Tumeur-ROI.nii']);

DCE_09  = MRIread([tumourPath,'M',num2str(rID),'_J09_Tumeur-DCE-ROI.nii.gz']);
DCE_11  = MRIread([tumourPath,'M',num2str(rID),'_J11_Tumeur-DCE-ROI.nii.gz']);
DCE_14  = MRIread([tumourPath,'M',num2str(rID),'_J14_Tumeur-DCE-ROI.nii.gz']);
DCE_16  = MRIread([tumourPath,'M',num2str(rID),'_J16_Tumeur-DCE-ROI.nii.gz']);


% 3) force time inclusion
dif = DCE_09.vol - DCE_11.vol;
DCE_11.vol(dif(:)>0) = 1;

dif = DCE_11.vol - DCE_14.vol;
DCE_14.vol(dif(:)>0) = 1;

dif = DCE_14.vol - DCE_16.vol;
DCE_16.vol(dif(:)>0) = 1;

% 4) force modalities inclusion
dif = DCE_09.vol - T2w_09.vol;
T2w_09.vol(dif(:)>0) = 1;

dif = DCE_11.vol - T2w_11.vol;
T2w_11.vol(dif(:)>0) = 1;

dif = DCE_14.vol - T2w_14.vol;
T2w_14.vol(dif(:)>0) = 1;

dif = DCE_16.vol - T2w_16.vol;
T2w_16.vol(dif(:)>0) = 1;

% 5) force time inclusion T2w
dif = T2w_09.vol - T2w_11.vol;
T2w_11.vol(dif(:)>0) = 1;

dif = T2w_11.vol - T2w_14.vol;
T2w_14.vol(dif(:)>0) = 1;

dif = T2w_14.vol - T2w_16.vol;
T2w_16.vol(dif(:)>0) = 1;


% % 6) save output
% MRIwrite(T2w_09,[tumourPath,'corrected/M',num2str(rID),'_J09_T2w.nii']);
% MRIwrite(T2w_11,[tumourPath,'corrected/M',num2str(rID),'_J11_T2w.nii']);
% MRIwrite(T2w_14,[tumourPath,'corrected/M',num2str(rID),'_J14_T2w.nii']);
% MRIwrite(T2w_16,[tumourPath,'corrected/M',num2str(rID),'_J16_T2w.nii']);
% 
% MRIwrite(DCE_09,[tumourPath,'corrected/M',num2str(rID),'_J09_DCE.nii']);
% MRIwrite(DCE_11,[tumourPath,'corrected/M',num2str(rID),'_J11_DCE.nii']);
% MRIwrite(DCE_14,[tumourPath,'corrected/M',num2str(rID),'_J14_DCE.nii']);
% MRIwrite(DCE_16,[tumourPath,'corrected/M',num2str(rID),'_J16_DCE.nii']);
% 


