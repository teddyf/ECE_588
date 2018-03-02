clear; close all

image_size = 100;
middle = floor(image_size/2);
N_images = 40;
radius = @(z) 6*z;
center = @(z) 3*[z+10,z];

B = make_dummy_image(image_size,radius,center,N_images);

figure(1); clf
for z=1:N_images
    imagesc(B(:,:,z));
    pause(0.2);
    colormap gray;
end