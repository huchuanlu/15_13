function [segNeib segNeibRat] = SegNeib(seg, nSeg)
% segNeib(i,j) = 1 if seg_i and seg_j are neighbors.
% segNeibRat(i,j) records the proportion of the boundary of the seg_j on the seg_i
nSeg = double(nSeg);
[dxl dxr dyt dyb] = Grad(seg);
ml = (dxl~=0); mr = (dxr~=0); mt = (dyt~=0); mb = (dyb~=0);
segNeib = zeros(nSeg);
for nSegIdx = 1:nSeg
    segMask = (seg==nSegIdx);

    l = dxl(segMask&ml)+nSegIdx; segNeib(nSegIdx,l) = 1;
    r = dxr(segMask&mr)+nSegIdx; segNeib(nSegIdx,r) = 1;
    t = dyt(segMask&mt)+nSegIdx; segNeib(nSegIdx,t) = 1;
    b = dyb(segMask&mb)+nSegIdx; segNeib(nSegIdx,b) = 1;

    dxl(segMask&ml) = l;
    dxr(segMask&mr) = r;
    dyt(segMask&mt) = t;
    dyb(segMask&mb) = b;
end

segNeibMap = max(max(dxl, dxr), max(dyt, dyb));
segNeibRat = zeros(nSeg);

for nSegIdx = 1:nSeg
    h = hist(segNeibMap(seg==nSegIdx), 0:double(nSeg));
    segNeibRat(nSegIdx,:) = h(2:end);
end
% shareLen = segNeibRat;
segNeibRat = diag(1./sum(segNeibRat,2))*segNeibRat;

return;


function [dxl dxr dyt dyb] = Grad(seg)

[h w] = size(seg);

dxl = zeros(h, w);   dxl(:, 2:w) = seg(:, 1:(w-1)) - seg(:, 2:w);
dxr = zeros(h, w);   dxr(:, 1:(w-1)) = seg(:, 2:w) - seg(:, 1:(w-1));
dyt = zeros(h, w);   dyt(2:h, :) = seg(1:(h-1), :) - seg(2:h, :);
dyb = zeros(h, w);   dyb(1:(h-1), :) = seg(2:h, :) - seg(1:(h-1), :);