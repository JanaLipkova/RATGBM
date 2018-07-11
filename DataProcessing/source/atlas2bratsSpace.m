%===========================================
%
% atlas 2 brats space
%
%============================================


% function atlas2bratsSpace

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')


%0) set path
dataPath = '../../atlas2brats/';

atlas = MRIread([dataPath,'atlas_pd.nii']);
brats = MRIread([dataPath,'t1_brats_space.nii']);

brats.vol = atlas.vol(end:-1:1, 1:end,1:end);
MRIwrite(brats,[dataPath,'atlas_pd_brats_space.nii']);
