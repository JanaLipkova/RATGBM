%===============================================
%
%   Convert mat data in given folder to nii
%
%   - data with ROI are already 3D matrices
%   - other data are 4D and need to be resize
%   - see ../doc/README for more info about the data
%
%  Copyright: Jana Lipkova
%  jana.lipkova@tum.de
%===============================================


% function convert_mat_anatomy_2nii

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')

% 1) Setup paths
rID = 42;
inputDataPath  = (['../../RAT_DATA_F98/M',num2str(rID),'/Anat-VOI/']);
outputDataPath = (['../../RAT_DATA_F98/M',num2str(rID),'/Anat-VOI/']);
inputNii       = ('../lib/T1w.nii');
FolderContent = dir(inputDataPath);


% 2) Read in some nifty to learn nifty structures
out = MRIread(inputNii);
data = out.vol;

% 3) Process data
% - read in files from the InputFolder
% - if ROI: save directetly to nii
% - else select correct data (no header), resize to 3D and save to nii
% - for T2 images, check if values above 1e+20 -> apply cliping

for i = 3:(size(FolderContent))
    
    % read file name from folder, create output file with nii extension, load data
    inputName   = FolderContent(i).name
    baseName    = inputName(1:end-4);
    outputName  = [outputDataPath, baseName,'.nii'];
    
    if(inputName(end-2:end) == 'mat')
        
        data        = load([inputDataPath,inputName]);
        bCollect = 0;
        
        
        dataOut = data.uvascroi().value;
        
        [Nx,Ny] =  size(dataOut);
        Nz = size(data.uvascroi,2);
        
        for iz = 1:Nz
            dataOut(:,:,iz) = data.uvascroi(iz).value;
        end;
        out.vol = dataOut;
        MRIwrite(out,outputName);
        
        
        %     clear data;
    end;
end
