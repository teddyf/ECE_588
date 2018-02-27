function images = make_dummy_image(image_size,radius,center,N_images)
images = zeros(image_size,image_size,N_images);
for z = 1:N_images
    radius_z = radius(z);
    center_z = center(z);
    images(:,:,z) = make_circle(image_size,radius_z,center_z);
end

end

function image = make_circle(image_size,radius,center)

image = zeros(image_size);
for i = 1:image_size
    for j = 1:image_size
        if((i-center(2)).^2+(j-center(1)).^2<radius^2)
            image(i,j) = 1;
        end
    end
end
end