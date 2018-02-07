%===================================
%
%   Crop Data
% -----------------------------------
%  Crop data by removing all non-zero elements based on the Segm mask
%
%  INPUT:
%     Segm = segmentation mask    
%     Vol  = medical volume
%
% OUTPUT:
%   SegmCrop = cropped segmentation mask
%   VolCrop  = cropped volume
%
%=================================== 

function VolCrop = cropData(Segm, Vol)

% 1st dim
mask = any(any(Segm));
Segm = Segm(:,:,mask);
Vol  = Vol(:,:,mask);

% Rotate to apply in 2nd dimension
Segm = rotate90_3D(Segm,2);
Vol = rotate90_3D(Vol,2);

mask = any(any(Segm));
Segm = Segm(:,:,mask);
Vol  = Vol(:,:,mask);


% Rotate to apply in 3rd dimension
Segm = rotate90_3D(Segm,3);
Vol = rotate90_3D(Vol,3);

mask = any(any(Segm));
Segm = Segm(:,:,mask);
Vol  = Vol(:,:,mask);

% Rotate to original coordinate system

Vol = rotate90_3D(Vol,2);
Vol = rotate90_3D(Vol,3);
Vol = rotate90_3D(Vol,3);
Vol = rotate90_3D(Vol,1);

Segm = rotate90_3D(Segm,2);
Segm = rotate90_3D(Segm,3);
Segm = rotate90_3D(Segm,3);
Segm = rotate90_3D(Segm,1);

% Vol = rotate90_3D(Vol,2);
% Segm = rotate90_3D(Segm,2);
% 
% Segm = permute(Segm,[2,1,3]);
% Vol = permute(Vol,[2,1,3]);

% Set return field
SegmCrop = Segm;
VolCrop = Vol;