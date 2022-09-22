function [ img ] = resizeHR(sigma,h,img,scale,method)
  img=imfilter(img,h,'circular');
  img=imresize(img,scale,method);
  noise=randn(size(img));
  img=img+sigma*noise;
end

