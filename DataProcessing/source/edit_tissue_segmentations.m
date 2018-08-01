% ===============================================
%  Binarise and Correct Tissue Segmentations
%
%
%
%  Copyright: Jana Lipkova
%  jana.lipkova@tum.de
% ===============================================

% function edit_tissue_segmentations

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')


%0) select case and set-up paths
rList = 42 %[29,30,34,36,38,42];


for rID = rList
    
    rID
    dataBase   = '../../RAT_DATA_F98/';
    
    %0) Set path to data
    anatomyPath = [dataBase,'M',num2str(rID),'/Anat-VOI/Registered/'];
    
    %1) Read-in mask tissue segmentaiotns
    mask =  MRIread([anatomyPath,'M',num2str(rID),'J09-Atlas-Mask_reg.nii']);
    wm   =  MRIread([anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_WM_reg.nii']);
    gm   =  MRIread([anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_GM_reg.nii']);
    csf   =  MRIread([anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_CSF_reg.nii']);
    
    % 2) Make mask binary
    mask.vol(mask.vol(:)>0.1) = 1;
    
    % 3) Make tissue binary, separate csf and tissue
    [Nx,Ny,Nz] = size(mask.vol);
    
    for iz=1:Nz
        for iy=1:Ny
            for ix=1:Nx
                
                if(mask.vol(ix,iy,iz) > 0)
                    
                    if(csf.vol(ix,iy,iz) > 0.4)
                        csf.vol(ix,iy,iz) = 1;
                        wm.vol(ix,iy,iz)  = 0;
                        gm.vol(ix,iy,iz)  = 0;
                    else
                        csf.vol(ix,iy,iz) = 0;
                        
                        if(wm.vol(ix,iy,iz) > 0.3)
                            wm.vol(ix,iy,iz) = 1;
                        end;
                        
                        if(gm.vol(ix,iy,iz) > 0.3)
                            gm.vol(ix,iy,iz) = 1;
                        end;
                        
                    end;
                    
                end;
            end;
        end;
    end;

    
    %4) Save ouptput
    MRIwrite(mask, [anatomyPath,'M',num2str(rID),'J09-Atlas_Mask_reg_cor.nii']);
    MRIwrite(gm,   [anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_GM_reg_cor.nii']);
    MRIwrite(wm,   [anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_WM_reg_cor.nii']);
    MRIwrite(csf,   [anatomyPath,'M',num2str(rID),'J09-VOI-Atlas_CSF_reg_cor.nii']);

end;