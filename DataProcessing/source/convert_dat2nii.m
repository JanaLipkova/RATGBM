%===================================
%
%   Convert folder with dat to nii
% -----------------------------------
%  INPUT:  path dat folder
%  OuTPUT: folder with nii
%
%===================================

% function convert_dat2nii

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/NIfTI_20140122/')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')


rID = 38;

% 1) Set up path
inputPath  = ['../../RAT_DATA_F98/M',num2str(rID,'%02d'),'/Results/MAP_Dc_10/'];
inputNii       = ('../lib/T1w.nii');

bResize = 1;


% 2) Read in some nifty to learn nifty structures
nii_data = MRIread(inputNii);


%4) i)  get names of files in the input folder
%   ii) convert them to nii
%   iii) save nii to output folder

files = dir(inputPath);
filesNames = {files.name};
filesNames = filesNames(~ismember(filesNames,{'.','..','.DS_Store'}));

for i = 1:length(filesNames)
    
    inFilename  = filesNames{i};
    outFilename = [inFilename(1:end-3),'nii'];
    
    if( inFilename(end-3:end) == '.dat')
        inFilename
        datVolume = loadMatrix([inputPath,inFilename]);
        
        
        if(bResize)
            newRes = [144,144,144];
            interp = '*linear';
            datVolume = resize3d(datVolume,  newRes,interp);
        end;
        
        nii_data.vol = datVolume;
        MRIwrite(nii_data, [inputPath,outFilename]);
    end
end;
