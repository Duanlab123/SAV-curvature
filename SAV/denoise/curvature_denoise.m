clear 
addpath util

sigma=20;
load lena.mat
uclean=double(im);
 
load lena_noise20.mat
[w,error,energy]=SAV_denoise(f);
psnr1=psnr(uint8(w),uint8(uclean)) 
figure;
plot(energy,'b','LineWidth',2);   
xlabel('Iteration')
ylabel('Energy')
 
figure;
imshow(uint8(w));