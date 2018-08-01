% ===============================================
%  Create brain mask
%
%   - used masked T2w scans to create brain masks
%
%
%  Copyright: Jana Lipkova
%  jana.lipkova@tum.de
% ===============================================

function threashold_gm_segmentations

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')


%0) select case and set-up paths
rList = [29,30,34,36,38,42];


for rID = rList
    
    rID
    dataBase   = '../../RAT_DATA_F98/';
    
    %0) Set path to data
    anatomyPath = [dataBase,'M',num2str(rID),'/Anat-VOI/Registered/'];
    
    %1) Read in tissue segmentaiotns
    gm  = MRIread([anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_GM_reg.nii']);
    MRIwrite(gm, [anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_GM_reg_orig.nii']);

    %3) Make the mask binary + save
    gm.vol(gm.vol(:)>0.2) = 1;
    
    %4) Save ouptput
    MRIwrite(gm, [anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_GM_reg.nii']);
    
end;