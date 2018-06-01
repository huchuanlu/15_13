function segment = SegmentImage(name)

sigma = 0.25;
kappa = 32;
nAppSP = 600;   % approximate number of superpixels

img = imread(name);
[h w nChannel] = size(img);
if(nChannel==1)
    img = repmat(img(:), [1 3]);
    img = reshape(img, h, w, 3);
end
nMinSize = int32(h*w/nAppSP);
[seg nSeg] = EGBIS(img, sigma, kappa, nMinSize);
segNeib = SegNeib(seg, nSeg);
[x, y] = meshgrid(1:nSeg, 1:nSeg);
mask = triu(segNeib,1);

segment.imgName = name;
segment.seg = seg;
segment.segArea = SegArea(seg, nSeg);
segment.pairSeg = int32([x(mask==1) y(mask==1)]);
segment.stat.P = nSeg;
segment.stat.E = int32(size(segment.pairSeg,1));