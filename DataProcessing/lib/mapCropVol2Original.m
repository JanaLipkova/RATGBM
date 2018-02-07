%===================================
%
%   Map cropped volume to Original
% -----------------------------------
%
%  INPUT:
%     orgBinary  =  original segmentation
%     cropBinary =  cropped segmentation   
%     cropVol    = the volume to be mapped based on the banry maps 
%
% OUTPUT:
%   outVol   = input volume mapped to original space based on input segmentations 
%
%=================================== 

function [outVol] = mapCropVol2Original(orgBinary,cropBinary,cropVol)

[Nx,Ny,Nz] = size(orgBinary);
[Cx,Cy,Cz] = size(cropBinary);

% 1st dim
mask = any(any(orgBinary));
tmpZ = find(mask(1,1,:)>0);
alphaZ = min(tmpZ);

% 2nd dim: Rotate to apply in 2nd dimension
orgBinary = rotate90_3D(orgBinary,2);

mask = any(any(orgBinary));
tmpY = find(mask(1,1,:)>0);
alphaY = Ny - max(tmpY) + 1;


% 3rd dim: Rotate to apply in 3rd dimension
orgBinary = rotate90_3D(orgBinary,3);

mask = any(any(orgBinary));
tmpX = find(mask(1,1,:)>0);
alphaX = min(tmpX);

outVol = zeros(Nx,Ny,Nz);
outVol(alphaX:alphaX+Cx-1, alphaY:alphaY+Cy-1,alphaZ:alphaZ+Cz-1) = cropVol;

 