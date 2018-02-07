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


% function convert_mat2nii

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')

% 1) Setup paths
inputDataPath  = ('../../Image_Analyses_data/');
outputDataPath = ('../../Image_Analyses_data_nii/');
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

for i = 698:(size(FolderContent)-1)
    i;
    
    % read file name from folder, create output file with nii extension, load data
    inputName   = FolderContent(i).name
    baseName    = inputName(1:end-4);
    outputName  = [outputDataPath, baseName,'.nii'];
    
    data        = load([inputDataPath,inputName]);
    bCollect = 0;
    
    if( baseName(end-2:end)=='ROI')
        bCollect=2;
    else
        if(baseName(8:10)=='VOI')
            bCollect=1;
        else
            if(baseName(9:13)=='Atlas')
                bCollect=2;
            end;
        end;
    end;
            
    if(bCollect==1)
        dataOut = data.uvascroi().value;
        
        [Nx,Ny] =  size(dataOut);
        Nz = size(data.uvascroi,2);
        
        for iz = 1:Nz
            dataOut(:,:,iz) = data.uvascroi(iz).value;
        end;
        out.vol = dataOut;
        MRIwrite(out,outputName);
    end;
    
    if (bCollect==2)
        id = str2num(baseName(2:3));
        if(id<29);
            Nz = 25;
        else
            Nz = 19;
        end;
        
        tmp = data.uvascroi().value;
        [Nx,Ny] =  size(tmp);
        dataOut = zeros(Nx,Ny,Nz);
        
        Nslices = size(data.uvascroi,2);
        
        for is = 1:Nslices
            sliceID = data.uvascroi(is).displayedslice;
            dataOut(:,:,sliceID) = data.uvascroi(is).value;
        end;
        out.vol = dataOut;
        MRIwrite(out,outputName);
        
    end;
    
     clear data;
end
