clear all
%%%%%%%%%%%%%  Data parameters
uclean=imread('img_001_mean.png');
f = imnoise(uclean,'poisson');
f= double(f)+5*randn(size(uclean));
[w,error,energy]=SAV_denoise(double(f));
figure;imshow(w,[0 255])
psnr1=psnr(uint8(w),uclean)


