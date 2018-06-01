clc;
clear all;
tic;
addpath('geObj_code');
addpath('seg_code');
%addpath('LLC');
on=3000;
%spnumber=100;% superpixel number
imgRoot='.\imgs\';% test image path
saldir=['.\out\'];% the output path of the saliency map
mkdir(saldir);
supdir='./superpixels/';% the superpixel label file path
nfeature=4;
if ~isdir(supdir)
mkdir(supdir);
end
if ~isdir(saldir)
mkdir(saldir);
end
imnames=dir([imgRoot '*' 'jpg']);
for ii=1:length(imnames)
    disp(ii);
    imname=[imgRoot imnames(ii).name]; 
    im=imread(imname);
    [w,h,dim]=size(im); %%the original size,
    if dim==1
        im0=zeros(w,h,3);
        im0(:,:,1)=im;
        im0(:,:,2)=im;
        im0(:,:,3)=im; 
        im=im0;
    end
%% cut the frame
[input_im,flag]=cutframe(im); 
input_im=uint8(input_im);
[m,n,k] = size(input_im);
%% objectness prior
[~, allboxes]=runObjectness(input_im,0);
[~,index]=sort(allboxes(:,5),'descend');
boxes=floor(allboxes(index(1:on),:));
mm0=zeros(m,n);
for i=1:on
   mm0(boxes(i,2):boxes(i,4),boxes(i,1):boxes(i,3))=mm0(boxes(i,2):boxes(i,4),boxes(i,1):boxes(i,3))+1;
end
if max(mm0(:))
mm0=normalize(mm0);
end
objectprior=saveimg(flag,mm0,[w,h],[saldir imnames(ii).name(1:end-4) '-object.png']);

imname=[imname(1:end-4) '.bmp'];
imwrite(input_im,imname);
smap=0;
%% save the codebook and the label
datap.rgb=[];
datan.rgb=[];
datap.lab=[];
datan.lab=[];
datap.lbp=[];
datan.lbp=[];
datap.hog=[];
datan.hog=[];
labelp=[];
labeln=[];
% save the test samples
rgbvals0=[];
labvals0=[];
hogvals0=[];
lbpvals0=[];
cutim0=zeros(m,n);
angle = 360;    bins = 8; 
[ vM iM ] = hog(rgb2gray(input_im), angle, bins);  
[A,~] = LBP_uniform(rgb2gray(input_im));
%% 
disp('building graph');
N = m*n;
imm = double((input_im));
m1=imm(:,:,1);
m2=imm(:,:,2);
m3=imm(:,:,3);
E = edges4connected(m,n);
V=1./(1+sqrt((m1(E(:,1))-m1(E(:,2))).^2+(m2(E(:,1))-m2(E(:,2))).^2+(m3(E(:,1))-m3(E(:,2))).^2));
AA=1000*sparse(E(:,1),E(:,2),0.3*V);
g = fspecial('gauss', [5 5], sqrt(5));

pos_mat=[];
pos_mat_b=[];
pos_mat_f=[];
%% 
count=0;
bsf0=[];
bsb0=[];
for ss=1  %% can be extended to multiscale
    temp_pos_mat=[];
    spnumber=(ss-0)*100;
%% ----------------------generate superpixels--------------------%%
   % the slic software support only the '.bmp' image
    comm=['SLICSuperpixelSegmentation' ' ' imname ' ' int2str(20) ' ' int2str(spnumber) ' ' supdir];
    system(comm);    
    spname=[supdir imnames(ii).name(1:end-4)  '.dat'];
    superpixels=ReadDAT([m,n],spname); % superpixel label matrix
    spnum=max(superpixels(:));% the actual superpixel number
% compute the feature (mean color in lab color space) 
% for each node (superpixels)
[m,n,k] = size(input_im);
    input_vals=reshape(input_im, m*n, k);
    rgb_vals=zeros(spnum,1,3);
    inds=cell(spnum,1);
    for i=1:spnum
        inds{i}=find(superpixels==i);
        rgb_vals(i,1,:)=mean(input_vals(inds{i},:),1);  
    end
%% color feature
lab_vals = colorspace('Lab<-', rgb_vals); 
lab_vals=reshape(lab_vals,spnum,3);
rgb_vals=reshape(rgb_vals,spnum,3);
 
     STA=regionprops(superpixels,'all');
       if ss==1
         STA0=STA;
       end
   lbp_vals=zeros(spnum,1,59);
 for i=1:spnum
        temp=A(STA(i).PixelIdxList);
        lbp_vals(i,1,:)=hist(temp,1:59);
 end
lbp_vals=reshape(lbp_vals,spnum,59);
%% hog feature
   hog_vals=zeros(spnum,bins);
  for i=1:spnum
 tempindex=STA(i).PixelIdxList;
  XX=iM(tempindex);
  YY=vM(tempindex);
  for ik=1:bins
      hog_vals(i,ik)=sum(YY(XX==ik));
  end
  if max(hog_vals(i,:))>0
  hog_vals(i,:)=normalize(hog_vals(i,:));
  end
  end
        
%%    
    % top, down, right, left
    bst=unique(superpixels(1,1:n));
    bsd=unique(superpixels(m,1:n));
    bsr=unique(superpixels(1:m,1));
    bsl=unique(superpixels(1:m,n));
    bs=sort(unique([bst,bsd,bsr',bsl'])); %index for image boundary superpixels
    spatial_prior=zeros(spnum,1);
    region_dist=zeros(spnum,spnum,4);
     [x y] = meshgrid(-n/2+1:n/2, -m/2+1:m/2); 
    x = x.^2;
    y = y.^2;
       k=2;
    centersp=setdiff([1:spnum],bs);
    for iii=1:spnum
        i=iii;
          temp_x_dist = mean(x(STA(i).PixelIdxList));
        temp_y_dist = mean(y(STA(i).PixelIdxList));
        temp_pos_mat =[temp_pos_mat; mean(x(STA(i).PixelIdxList)), mean(y(STA(i).PixelIdxList))]; 
        spatial_prior(i) = exp(-k^2*temp_x_dist/n^2 - k^2*temp_y_dist/m^2); %spatial weight
        for ix=1:length(bs)
            region_dist(i,bs(ix),1)=histDist(hog_vals(i,:),hog_vals(bs(ix),:));
            region_dist(i,bs(ix),2)=sum((rgb_vals(i,:)-rgb_vals(bs(ix),:)).^2);
            region_dist(i,bs(ix),3)=sum((lab_vals(i,:)-lab_vals(bs(ix),:)).^2);
            region_dist(i,bs(ix),4)=histDist(lbp_vals(i,:),lbp_vals(bs(ix),:));
        end
    end
    for i=1:4  %normalization within the four feature spaces
    region_dist(:,:,i)=(region_dist(:,:,i)-min(min(region_dist(:,:,i))))/(max(max(region_dist(:,:,i)))-min(min(region_dist(:,:,i))));
    end

   region_dist= sum(region_dist,3);   
  region_dist=region_dist+region_dist'; 
  temp_sal = abs(-log(1 - region_dist+0.0000001));%avoid the complex number of the log computation

  objectP=zeros(spnum,1); %objectprior
  for i=1:spnum
     objectP(i,1)= mean(mm0(STA(i).PixelIdxList)); 
  end
   darkpa00=objectP;
   objectP=exp(objectP);
 saliency = spatial_prior .* sum(temp_sal(:,bs), 2).*objectP;
 map=zeros(m,n);
   for i=1:spnum
      map(STA(i).PixelIdxList)=repmat(saliency(i), size(STA(i).PixelIdxList)); 
   end
   map(map==inf)=max(map(:));
   %%
 %% graphcut 
[cutim]=graphcut0(AA,g,map); 
 %cutim0=normalize(cutim0+cutim);
 maxi=0;
 cutim0=cutim+cutim0;
 %% MKB
    bsf=[];
    bsb=[];
    meancut=1.5.*mean(cutim(:));
    thresh=0.05;
  %  maxi=0;
    for i=1:spnum 
        meantemp=mean(cutim(STA(i).PixelIdxList));
        if meantemp>maxi
            maxi=meantemp;
            maxind=i;
        end
        if meantemp>= meancut
            bsf=[bsf,i];
        else if meantemp< thresh
            bsb=[bsb,i];
            end
        end 
    end
    if isempty(bsf)
       bsf=[maxind]; 
    end
       bsb=unique([bsb,bs]);
   
 datap.rgb=[datap.rgb;rgb_vals(bsf,:)];
 datan.rgb=[datan.rgb;rgb_vals(bsb,:)];
 datan.hog=[datan.hog;hog_vals(bsb,:)];
 datap.hog=[datap.hog;hog_vals(bsf,:)];
 datap.lab=[datap.lab;lab_vals(bsf,:)];
 datan.lab=[datan.lab;lab_vals(bsb,:)];
 datap.lbp=[datap.lbp;lbp_vals(bsf,:)];
 datan.lbp=[datan.lbp;lbp_vals(bsb,:)];
 
labelp=[labelp;repmat(1,length(bsf),1)];
labeln=[labeln;repmat(-1,length(bsb),1)];
rgbvals0=[rgbvals0;rgb_vals];
labvals0=[labvals0;lab_vals];
hogvals0=[hogvals0;hog_vals];
lbpvals0=[lbpvals0;lbp_vals];


pos_mat=[pos_mat;temp_pos_mat];
pos_mat_b=[pos_mat_b;temp_pos_mat(bsb,:)];
pos_mat_f=[pos_mat_f;temp_pos_mat(bsf,:)];

bsb0=[bsb0,bsb+count];
bsf0=[bsf0,bsf+count];
count=spnum+count;

end

data.rgb=[datap.rgb;datan.rgb];
data.lab=[datap.lab;datan.lab];
data.hog=[datap.hog;datan.hog];
data.lbp=[datap.lbp;datan.lbp];
 
label.rgb = [labelp;labeln];
label.lab = label.rgb;
label.hog = label.rgb;
label.lbp = label.rgb;

datanum = [size(data.rgb,1) size(data.lab,1)];
data0.rgb=rgbvals0;
data0.lab=labvals0;
data0.hog=hogvals0;
data0.lbp=lbpvals0;

%% LLC
beta=0.1;
knn=floor(size(datan.rgb,1)*2/3);

datan.rgb=[datan.rgb,datan.lab,datan.hog,datan.lbp];
data0.rgb=[data0.rgb,data0.lab,data0.hog,data0.lbp];
Coeff=LLC_coding_dis1(datan.rgb,data0.rgb,pos_mat,pos_mat_b, knn, beta);
ReError=sum((data0.rgb-Coeff*datan.rgb).^2,2);
Rec0=(ReError- min(ReError(:)))/(max(ReError(:)) - min(ReError(:)));

knn=floor(size(datap.rgb,1)*2/3);

datap.rgb=[datap.rgb,datap.lab,datap.hog,datap.lbp];
Coeff=LLC_coding_dis1(datap.rgb,data0.rgb,pos_mat,pos_mat_f, knn, beta);
ReError=sum((data0.rgb-Coeff*datap.rgb).^2,2);
Rec1=1-(ReError- min(ReError(:)))/(max(ReError(:)) - min(ReError(:)));

Rec=(Rec0.*Rec1);
map0=zeros(m,n);  
spnum0=size(STA0,1);

 conf=Rec;
  for i=1:spnum0
      map0(STA0(i).PixelIdxList)=repmat(conf(i), size(STA0(i).PixelIdxList));  
  end
map0=normalize(map0);
[mkbcut0]=graphcut0(AA,g,map0); 
mkbcut=normalize(mkbcut0);
step2graph=saveimg(flag,mkbcut,[w,h],[saldir imnames(ii).name(1:end-4) '-step2.png']);
%% fusion of the local and global map
cutim=normalize(cutim0);
cutim=saveimg(flag,cutim,[w,h],[saldir imnames(ii).name(1:end-4) '-step1.png']);
sal=normalize(cutim+step2graph);
imwrite(sal,[saldir imnames(ii).name(1:end-4)  '-final.png']); 
end







