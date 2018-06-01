function maskImg = SegImg_display(segment, disObjLine)
if(~exist('disObjLine', 'var'))
    disObjLine = 'r';
end
img = imread(segment.imgName);
seg = segment.seg;

maskImg = img;
maskImgR = img(:,:,1);
maskImgG = img(:,:,2);
maskImgB = img(:,:,3);

[dxl, ~, dyt] = Grad(seg);
if(strcmpi(disObjLine, 'r'))
    maskImgR((dxl~=0)|(dyt~=0)) = 255;
    maskImgG((dxl~=0)|(dyt~=0)) = 0;
    maskImgB((dxl~=0)|(dyt~=0)) = 0;
elseif(strcmpi(disObjLine, 'b'))
    maskImgR((dxl~=0)|(dyt~=0)) = 0;
    maskImgG((dxl~=0)|(dyt~=0)) = 0;
    maskImgB((dxl~=0)|(dyt~=0)) = 255;
end

maskImg(:,:,1) = maskImgR;
maskImg(:,:,2) = maskImgG;
maskImg(:,:,3) = maskImgB;
return;

function [dxl dxr dyt dyb] = Grad(seg)

[h w] = size(seg);

dxl = zeros(h, w);   dxl(:, 2:w) = seg(:, 1:(w-1)) - seg(:, 2:w);
dxr = zeros(h, w);   dxr(:, 1:(w-1)) = seg(:, 2:w) - seg(:, 1:(w-1));
dyt = zeros(h, w);   dyt(2:h, :) = seg(1:(h-1), :) - seg(2:h, :);
dyb = zeros(h, w);   dyt(1:(h-1), :) = seg(2:h, :) - seg(1:(h-1), :);