function [ du ] = downsa( u,scale,kernel )
du=imresize(u, 1/scale,kernel);
end

