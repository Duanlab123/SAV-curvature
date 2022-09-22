clear all 
addpath util
%kernel = fspecial('motion',20,pi/4);
kernel = fspecial('gaussian',[5 5],5);
uclean=double(imread('kodim14.png'));
result=zeros(size(uclean));
 for Nc=1:3
   img=uclean(:,:,Nc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   f=imfilter(img,kernel,'circular','conv');  
   f=f+5*randn(size(f));
   [w,error,Energy_iter]=SAV_deblur(f,kernel);
   result(:,:,Nc)=w;
 end
psnr1=psnr(uint8(result),uint8(uclean))
figure;imshow(uint8(result),[])


