clear all;clc;close all;
% % line 498 598 718
datamount = 245760;
tic
time=datetime('now','Format','yMd_HH,mm,ss');
timec=char(time);
mkdir ('./fftplate',timec)
mkdir ('./output')
load('nyu_depth_v2_labeled.mat', 'depths')
load('nyu_depth_v2_labeled.mat', 'images')

filename = 'val_list.txt'
% normalize = true;
filepath = './result';
midname = '/result_';

filepath2 = './paper_results';
p1 = '/CSPN_results/eval_result';
p2 = '/NLSPN_NYU_results';
filepath3 = './examples';
midname2 = '/rgb_';
idx = readmatrix(filename);
% I_our = {};
I_hed = {};
I_rgb = {};
our_missed= [];
I_depth = {};
for i = 1:length(idx)
    if isfile([filepath midname num2str(i) '.png'])

        I_hed{i} = imread([filepath3 midname2 num2str(idx(i)) '.png']);
        I_rgb{i} = images(:,:,:,idx(i));
        I_depth{i} = 1000 * double (depths(:,:,idx(i))); 
    else

         I_hed{i}=[]; 
        our_missed = [our_missed;i];
    end
   
end











for sss = 1:length(I_depth) 
if length(I_depth{sss}) ~= 0

% % % % % % % % % % % 
I = I_hed{sss};
I=rgb2gray(I);


depthimage = I_depth{sss};
rgb=I_rgb{sss};
% % 698 line, random 1500 points
ranum = randi([1 length(reshape(depthimage,[],1))],1,datamount);
[row_dd col_dd] =  ind2sub(size(depthimage),ranum);
depthimage_temp = zeros(size(depthimage));
for ddd = 1:length(row_dd)

depthimage_temp(row_dd(ddd),col_dd(ddd)) = depthimage(row_dd(ddd), col_dd(ddd));
end
depthimage = depthimage_temp;

level=graythresh(I);

 BW = imbinarize(I,level);






[B,L] = bwboundaries(BW,'holes');
fg=label2rgb(L, @jet, [.5 .5 0]);

%depth image to point cloud
tptcloud=depthToCloud(depthimage,[0 0 0]);%organize point cloud
ptcloud=pointCloud(tptcloud,'color',rgb);
% ptcloud = pcdenoise(ptcloud);
tptcloudr=reshape(tptcloud,[],3);%reshape 
ptcloudr=pointCloud(tptcloudr);

depthimage2=depthimage;
tptcloud2=depthToCloud(depthimage2,[0 0 0]);
pc3=pointCloud(tptcloud2,'color',rgb);
% figure
% pcshow(pc3)
depthimage=depthimage2;


%depth image to point cloud
tptcloud=depthToCloud(depthimage,[0 0 0]);
%organize point cloud
ptcloud=pointCloud(tptcloud,'color',rgb);

tptcloudr=reshape(tptcloud,[],3);%reshape 
ptcloudr=pointCloud(tptcloudr);
%
rgbr=reshape(rgb,[],3);

degree=8;

 A = [1 0 0 0; ...
     0 1 0 0; ...
     0 0 1 0; ...
     0 0 0 1];
tform = affine3d(A);
ptcloudrot = pctransform(ptcloud,tform);

%raw parameter
xmin=-0.3;
xmax=0.5;
ymin=-0.06;
ymax=0.44;
zmin=1.7;
zmax=2.8;

pttemp=ptcloudrot.Location;
 pttempc=[reshape(pttemp,[],3) double(rgbr)];
%delete outlier  (min to max )
pttempc=[reshape(pttemp,[],3) double(rgbr)];

sizep=size(pttemp);
pttp=reshape(pttempc(:,1:3),sizep(1,1),sizep(1,2),sizep(1,3));
pttc=reshape(pttempc(:,4:6),sizep(1,1),sizep(1,2),sizep(1,3));
% tpt=pttempc(:,1:3);
% tc=uint8(pttempc(:,4:6));
ptcleartable=pointCloud(pttp,'color',uint8(pttc));


[temp1c,temp1r]=find(isnan(pttp(:,:,1)));%Because is depth image, once x is nan, y and z bound to be nan.





%///////////////////////////////////////////

 cform = makecform('srgb2lab');
lab_he = applycform(he,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

nColors = length(B);
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);
pixel_labels = reshape(cluster_idx,nrows,ncols);

segmented_images = cell(1,3);% squeeze the cluster
rgb_label = repmat(pixel_labels,[1 1 3]);

for k = 1:nColors
    color = he;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
end
% find accurate red!!
fi1=segmented_images{1};
fi2=segmented_images{2};
fi3=segmented_images{3};
fi1ri=find(fi1(:,:,1)>230);
fi2ri=find(fi2(:,:,1)>230);
fi3ri=find(fi3(:,:,1)>230);
length1=length(fi1ri);
length2=length(fi2ri);
length3=length(fi3ri);
for i=1:3
    if(length1>length2 && length1>length3)
       redc=fi1;
    elseif(length2>length1 && length2>length3)
        redc=fi2;
    else
        redc=fi3;
    end
end
% redc is ture red
%///////////////////////////////////////
% figure
% imshow(redc);title('objects in cluster 2')

for i=1:length(segmented_images)
    BWALL1=rgb2gray(segmented_images{i});
        BWALL2=segmented_images{i};
% %      removes the background and edge image
if(BWALL2(:,:,1)==BWALL2(:,:,2) )
    continue
end
BWALL{i}=imbinarize(BWALL1);
end
colorfinal(1:sizep(1),1:sizep(2),1)=0;
colorfinal(1:sizep(1),1:sizep(2),2)=0;
colorfinal(1:sizep(1),1:sizep(2),3)=0;
cont=0;
BWALLf=cell(1);

exist_rate_set=0.52;
raw_rate=[];
sizecheck=size(depthimage);

for yai=1:length(BWALL)
    critzero=zeros(sizecheck);
    crit=BWALL{yai};
    critre=reshape(crit,[],1);
    for1=find(critre==1);
    critzero(for1)=depthimage(for1);
    inifor1re=reshape(logical(critzero),[],1);
    depfor1=find(inifor1re==1);
    exist_rate=(length(depfor1)/length(for1));
    raw_rate=[raw_rate;exist_rate];
    if (isnan(exist_rate) || exist_rate==0) 
        continue
    end
    if (exist_rate>=exist_rate_set)
        cont=cont+1;
        BWALLf{cont}=BWALL{yai};
        continue
    else
        indicemount_for1=length(for1);
        for2=for1(1:floor((indicemount_for1)/2),1);
        for3=for1(floor((indicemount_for1)/2)+1:length(for1)-1,1);
    end
    % for2
    critzero(for2)=depthimage(for2);
    inifor1re=reshape(logical(critzero),[],1);
    depfor1=find(inifor1re==1);
    exist_rate2=(length(depfor1)/length(for2));
    
    
    if(exist_rate2>exist_rate_set)
        cont=cont+1;
        critzero=zeros(sizecheck);
        critzero(for2)=crit(for2);
        BWALLf{cont}=critzero;
        
    else
        indicemount_for2=length(for2);
        for4=for2(1:floor((indicemount_for2)/2),1);
        for6=for2(floor((indicemount_for2)/2)+1:length(for2)-1,1);
    end
    %3
         critzero(for3)=depthimage(for3);
    inifor1re=reshape(logical(critzero),[],1);
    depfor1=find(inifor1re==1);
    exist_rate3=(length(depfor1)/length(for3));
    
    
    if(exist_rate3>exist_rate_set)
        cont=cont+1;
                critzero=zeros(sizecheck);
        critzero(for3)=crit(for3);
        BWALLf{cont}=critzero;
        continue
    else
        indicemount_for3=length(for3);
        for5=for3(1:floor((indicemount_for3)/2),1);
        for7=for3(floor((indicemount_for3)/2)+1:length(for3)-1,1);
    end
        critzero(for4)=depthimage(for4);
        inifor1re=reshape(logical(critzero),[],1);
        depfor1=find(inifor1re==1);
        exist_rate4=(length(depfor1)/length(for4));
        if (exist_rate4>exist_rate_set)
             cont=cont+1;
                     critzero=zeros(sizecheck);
        critzero(for4)=crit(for4);
        BWALLf{cont}=critzero;
        end
        critzero(for6)=depthimage(for6);
        inifor1re=reshape(logical(critzero),[],1);
        depfor1=find(inifor1re==1);
        exist_rate6=(length(depfor1)/length(for6));
        if (exist_rate6>exist_rate_set)
             cont=cont+1;
                     critzero=zeros(sizecheck);
        critzero(for6)=crit(for6);
        BWALLf{cont}=critzero;
        end
    %3t

    
        critzero(for5)=depthimage(for5);
        inifor1re=reshape(logical(critzero),[],1);
        depfor1=find(inifor1re==1);
        exist_rate5=(length(depfor1)/length(for5));
        if (exist_rate5>exist_rate_set)
             cont=cont+1;
                     critzero=zeros(sizecheck);
        critzero(for5)=crit(for5);
        BWALLf{cont}=critzero;
        end
        critzero(for7)=depthimage(for7);
        inifor1re=reshape(logical(critzero),[],1);
        depfor1=find(inifor1re==1);
        exist_rate7=(length(depfor1)/length(for7));
        if (exist_rate7>exist_rate_set)
             cont=cont+1;
                     critzero=zeros(sizecheck);
        critzero(for7)=crit(for7);
        BWALLf{cont}=critzero;
        end
end
        


        


        




% % % % Running!
for yay=1:length(BWALLf)
    redc=uint8(BWALLf{yay});
    redc_size=size(redc);
    if (redc_size(1)~=sizecheck(1) && redc_size(2)~=sizecheck(2))
        continue
    end
rgbhal=redc;
[ar ac]=find(redc(:,:,1)>0);
% loading region


t=[ar ac];
stop_critiria=length(t);
if (length(t)==0)
    continue
end
check1=sub2ind(size(depthimage2),t(:,1),t(:,2));
checkpoint1=find(depthimage(check1)==0);

if(stop_critiria<3500 || length(t)==length(checkpoint1))
    continue
    
end


%indice of region
psumx=0;
psumy=0;
psumz=0;
o=0;
pttpnanin=[];
pttind=[];
pttindgre=[];
for k=1:length(t)
    redc(t(k,1),t(k,2),1)=80;
    if(~isnan(pttp(t(k,1),t(k,2),1)))
        o=o+1;
    ptemx=pttp(t(k,1),t(k,2),1);
    ptemy=pttp(t(k,1),t(k,2),2);
    ptemz=pttp(t(k,1),t(k,2),3);
    psumx=psumx+ptemx;
    psumy=psumy+ptemy;
    psumz=psumz+ptemz;
    pttindgre=[pttindgre;t(k,1) t(k,2)];
    else
        % kinect 1 pixel to real distance is 6*4.5mm, f=2.9mm(rgb) f=6.1mm(IR)
        tx=t(k,1)/sizep(1,1)*(4.5*10^-4);
         ty=t(k,1)/sizep(1,2)*(6*10^-4);
        pttpnanin=[pttpnanin;tx ty  ];
        pttind=[pttind;t(k,1) t(k,2)];
    end
    rgbrealcir(t(k,1),t(k,2),:)=rgb(t(k,1),t(k,2),:);
end
%display fix point on image
for k=1:length(pttind)
    redc(pttind(k,1),pttind(k,2),3)=200;
    redc(pttind(k,1),pttind(k,2),1)=49;
end

A=[];
B=[];
for o=1:length(pttind)
a1=ismember(pttindgre(:,1),pttind(o,1));
b1=ismember(pttindgre(:,2),pttind(o,2));
ij=a1+b1;
iji=find(ij==2);
pttindgre2=pttindgre;
pttindgre(iji,:)=[];
end



for i=1:length(pttindgre)
    un(pttindgre(i,1),pttindgre(i,2))=depthimage(pttindgre(i,1),pttindgre(i,2));

end

sizeun=size(un);

un2(sizeun(1,1)+1:2*sizeun(1,1),sizeun(1,2)+1:2*sizeun(1,2))=0;
un2(1:sizeun(1,1),1:sizeun(1,2))=un;
un3=un2;%store un2 

ori=un2;


cycle=1500;% loop cycle
loop=0;
inini=(pttindgre(:,2)-1)*480+pttindgre(:,1);
inini=find(un2~=0);

sizeun2=size(un2);
rx=max(t(:,2))-min(t(:,2));
ry=max(t(:,1))-min(t(:,1));

rinitial=[sizeun2(1,1) sizeun2(1,2)]./[ry rx];
% % extract spatial frequency from rgb 2D data
% gray level first......


rgbg=rgb2gray(rgbrealcir);

rgb2=imresize(rgbg,[sizeun2(1,1) sizeun2(1,2)]);
frergb=fft2(double(rgb2));
F2 = fftshift(frergb);
frergbabs1=abs(F2);

frergbabsDISPLAY = log(frergbabs1+1);
frergbabsDISPLAY= mat2gray(frergbabsDISPLAY); % Use mat2gray to scale the image between 0 and 1
fg=ones(sizeun2(1),sizeun2(2));


% %  Morphologically close method

f2d=(abs(frergbabsDISPLAY));
f2d2=double(f2d);
f2d2=reshape(f2d2,[],1);
mrgb=find(f2d2==max(f2d2));% % m is the 0 frequency of picture.
f2d2(mrgb)=0;
n=find(f2d2==max(f2d2));% % n is the main order of fft

se = strel('disk',3,0);% 3rd number to be zero to be tight to white area.

F2r=reshape(f2d,[],1);
sizef2d=size(f2d);
f2d=mapminmax(F2r,0,1);
f2d=reshape(f2d,sizef2d(1),sizef2d(2));

IM2=imdilate(f2d,se);


sizeIM2=size(IM2);
IM2r2=reshape(IM2,[],1);
IM2r2=abs(IM2r2);
IM2r2=mapminmax(IM2r2,0,1);
IM2=reshape(IM2r2,sizeIM2(1),sizeIM2(2));

frergbabs=IM2;
fre2=frergbabsDISPLAY;
fre3=fre2;
fre2=double(fre2);
fre2=reshape(fre2,[],1);
m=find(fre2==max(fre2));% % m is the 0 frequency of picture.
fre2(m)=0;
n=find(fre2==max(fre2));% % n is the main order of fft
fre2(n)=0;
o=find(fre2==max(fre2));
[rowm colm] = ind2sub(size(fre3), m);
[rown coln] = ind2sub(size(fre3), n);
[rowo colo] = ind2sub(size(fre3), o);
% % cut 0 order

frergbabs=abs(frergbabs);
level2=graythresh(frergbabs);


rc=max(reshape([rowm-rown colm-coln],[],1));



% a(n)=255;
sizea=size(fre3);
fre2=reshape(fre2,sizea(1,1),sizea(1,2));





% % Display ft plane


ksr=struct([]);%clustering
g=1;% for ksr struct
ksr{1}=1;

ini=1;
% frergbabs=frergbabs/255;
for k=1:cycle
% % display ft plane
hj=find(isnan(un2));
un2(hj)=0;
fre=fft2(un2);
F = fftshift(fre); % Center FFT
fftresult=F;

%Circle with radius centered at 40,40 in an 80x80image

a=sizeun2(1,1)-1;
b=sizeun2(1,2)-1;
% r2= 5.0241+0.04*(k-1);% RADIUS
% r2=rinitial+(sizeun2(1,2)/cycle)*(k-1);% RADIUS

rb=0.5*rinitial(1,2)+(sizeun2(1,2)/(15*cycle))*(k-1);% RADIUS . rb in x, ra in y; 30 to be the best Trial
ra=0.5*rinitial(1,1)+(sizeun2(1,1)/(15*cycle))*(k-1);

[x y]=meshgrid(-b/2:b/2,-a/2:a/2);
C=zeros(a+1,b+1);

C((x/rb).^2+(y/ra).^2<=1)=1;

if ( k<100)
filitered=fftresult.*frergbabs.*C;
else
    filitered=fftresult.*C;
end

filitereddisplay=filitered;

filitered=ifftshift(filitered);
filiterresult=ifft2(filitered);


un2=abs(filiterresult);
un4=un2;
% %  store fft plate
    if (mod(k,100)==0)
        filitereddisplay=abs(filitereddisplay);
        filitereddisplay = log(filitereddisplay+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
        filitereddisplay= mat2gray(filitereddisplay); % Use mat2gray to scale the image between 0 and 1
     
        ksr{ini}=filitereddisplay;
         filenametemp=['fftplate/',timec,'/[ ',num2str(yay),' ]%d.png'];
         fil=char(filenametemp);
        fileName = sprintf(fil,k);
        
        imwrite (ksr{ini}, fileName);
        ini=ini+1;
%         rgbImage = ind2rgb(grayImage, jet(256));
    end

un2(inini)=un3(inini);
% end
    
end
absresult=abs(filiterresult);

un2(inini)=un3(inini);

un5=un3;

for p=1:length(pttind)
  depthimage(pttind(p,1),pttind(p,2))=un2(pttind(p,1),pttind(p,2));
  colorfinal(pttind(p,1),pttind(p,2),3)=25;
  colorfinal(pttind(p,1),pttind(p,2),1)=212;
  colorfinal(pttind(p,1),pttind(p,2),2)=212;
end
for p=1:length(pttind)
  depthimage(pttind(p,1),pttind(p,2))=un2(pttind(p,1),pttind(p,2));
end
% clear un2 un un3 C
end
F_abs=abs(F);
F_abs = log(F_abs+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
F_abs= mat2gray(F_abs); % Use mat2gray to scale the image between 0 and 1

filenametemp2=['fftplate/',timec,'/[ ',num2str(yay),' ]result.mat'];
         fil2=char(filenametemp2);
% depthFFT=un2;
morphology_mask=frergbabs;

% save(fil2,'depthimage','morphology_mask','ksr','colorfinal')
end
imwrite(uint16(depthimage),['./output/result_' num2str(sss) '.png'])
clear un2 un un3 C
end

toc
