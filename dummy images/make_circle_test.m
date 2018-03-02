clear; close all

image_size = 1000;
middle = floor(image_size/2);
N_images = 100;
radius = @(z) z;
center = @(z) 0.5*[z+5,z];

B = make_dummy_image(image_size,radius,center,N_images);

figure(1); clf
for z=1:N_images
    imagesc(B(:,:,z));
    pause(0.2);
    colormap gray;
end

STLname = 'dummySTL';
gridDATA = B;
STLformat = 'binary';
gridX = 1:image_size;
gridY = 1:image_size;
gridZ = 1:N_images;
CONVERT_voxels_to_stl(STLname,gridDATA,gridX, gridY, gridZ,STLformat);