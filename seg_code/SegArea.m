function segArea = SegArea(seg, nSeg)
segArea = zeros(nSeg,1);
seg = seg(:);
for k=1:nSeg
    segArea(k) = sum(seg==k);
end