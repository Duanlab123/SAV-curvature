clear all 
addpath util
uclean=double(imread('kodim23.png'));
upscaling=4;
kernel='bicubic';
load('h')
HR=zeros(size(uclean));
for Nc=1:3
  img=uclean(:,:,Nc);
  img=modcrop(img,upscaling);
  sigma=5;%nois in SR_mode
  low=resizeHR(sigma,h,img,1/upscaling,kernel);    
  HR(:,:,Nc)=SAV_SR(low,upscaling,kernel,h);
  wang_p=psnr(uint8(img),uint8( HR(:,:,Nc)))
  figure;imshow(HR(:,:,Nc),[])
 end
psnr1=psnr(uint8(HR),uint8(uclean))
figure;imshow(uint8(HR),[])
    

 