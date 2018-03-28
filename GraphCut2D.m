% tunable values
start_slice = 200;
transparency_level = 0.6;
pixel_delta = 2;
foreground_sampling = 0.2;
background_sampling = 0.2;
mask_depth = 50;

[imaVOL,scaninfo] = loadminc('group2/06/06_mr_tal.mnc');

vol = zeros(size(imaVOL));

% convert all slices to gray images
for slice = 1:size(imaVOL,3)
    vol(:,:,slice) = mat2gray(int16(imaVOL(:,:,slice)));
end

% segment out a single slice
slice_mid = vol(:,:,start_slice);
RGB = slice_mid;

% have user manually choose foreground and background points
imshow(RGB);
[f_x,f_y] = ginput;
[b_x,b_y] = ginput;

% ensure these are rounded (note: this was to fix a weird error I
% was getting, but not sure if it is actually necessary)
f_x = ceil(f_x);
f_y = ceil(f_y);
b_x = ceil(b_x);
b_y = ceil(b_y);

% convert pixel coords to indices
f_indices = sub2ind(size(RGB),f_y,f_x);
b_indices = sub2ind(size(RGB),b_y,b_x);
Mask = calcMask(RGB, f_indices, b_indices);
AlphaLevels = (1-Mask)*(1-transparency_level)+transparency_level;

% mark points
RGB_marked = markRGBwithFandBPoints(RGB, f_indices, b_indices);

% show these results
figure;
h = imshow(RGB_marked);
set(h, 'AlphaData', AlphaLevels);

vol_masked = vol;

source_Mask = Mask;
Mask_all = zeros(size(imaVOL));
for delta = -1:-1:-mask_depth
    start_slice+delta
    % calculate things for next slice
    slice_p = vol(:,:,start_slice+delta);
    RGB_p = slice_p;
    [f_indices_p, b_indices_p, RGB_p] = calcNewIndices(Mask, RGB_p, pixel_delta, foreground_sampling, background_sampling);
    
    if size(f_indices_p,1) ~=0
        Mask_p = calcMask(RGB_p, f_indices_p, b_indices_p);
        AlphaLevels_p = (1-Mask_p)*(1-transparency_level)+transparency_level;

        % now mark used coordinates
        RGB_p_marked = markRGBwithFandBPoints(RGB_p, f_indices_p, b_indices_p);

        % show figure
        %figure;
        %h = imshow(RGB_p_marked);
        %set(h, 'AlphaData', AlphaLevels_p);

        vol_masked(:,:,start_slice+delta) = slice_p + Mask_p;
        %+ (Mask_p.*0.5);

        Mask = Mask_p;
        Mask_all(:,:,start_slice+delta)= Mask_p;
    end
end
Mask = source_Mask;

% evolve up from source
for delta = 1:mask_depth
    start_slice+delta
    % calculate things for next slice
    slice_p = vol(:,:,start_slice+delta);
    RGB_p = slice_p;
    [f_indices_p, b_indices_p, RGB_p] = calcNewIndices(Mask, RGB_p, pixel_delta, foreground_sampling, background_sampling);
    
    if size(f_indices_p,1) ~=0
        Mask_p = calcMask(RGB_p, f_indices_p, b_indices_p);
        AlphaLevels_p = (1-Mask_p)*(1-transparency_level)+transparency_level;

        % now mark used coordinates
        RGB_p_marked = markRGBwithFandBPoints(RGB_p, f_indices_p, b_indices_p);

        % show figure
        figure;
        h = imshow(RGB_p_marked);
        set(h, 'AlphaData', AlphaLevels_p);

        vol_masked(:,:,start_slice+delta) = slice_p + Mask_p;
        %+ (Mask_p.*0.5);

        Mask = Mask_p;
    end
end

%% Crop 3D mask to tumor location

N_slices = size(vol_masked,3);
figure(1); clf
Mask_all = logical(Mask_all);
count = 0;
NZ_slices = [];

for i = 1:N_slices
    slice = Mask_all(:,:,i);
    if(sum(slice(:))~=0)
        NZ_slices = [NZ_slices,i]; %keep non-zero slices of mask
    end
end
cropped_im = Mask_all(:,:,NZ_slices); %cropped mask

%% Make STL file
[N1,N2,N3] = size(cropped_im);
STLname = '50slicesSTL.stl';
gridDATA = cropped_im;
STLformat = 'binary';
gridX = 1:N1;
gridY = 1:N2;
gridZ = 1:N3;
CONVERT_voxels_to_stl(STLname,gridDATA,gridX, gridY, gridZ,STLformat);


function [f_indices_p, b_indices_p, RGB_p] = calcNewIndices(Mask, RGB, pixel_delta,f_sampling, b_sampling)

    f_x_2 = [];
    f_y_2 = [];
    b_x_2 = [];
    b_y_2 = [];
    
    [rows,cols,dim] = size(Mask);
    
    for row = 1:rows
        for col = 1:cols
            element = Mask(row,col);
            X = rand;
            if element == 1
                if X < f_sampling
                    if pixelMeetsRequirements(pixel_delta,Mask,row,col)
                        f_x_2 = [f_x_2,row];
                        f_y_2 = [f_y_2,col];
                    end
                end
            end
            if element == 0
                if X < b_sampling
                    if pixelMeetsRequirements(pixel_delta,Mask,row,col)
                        b_x_2 = [b_x_2,row];
                        b_y_2 = [b_y_2,col];
                    end
                end
            end
        end
    end

    RGB_p = RGB;
    f_indices_p = transpose(sub2ind(size(RGB),f_x_2,f_y_2));
    b_indices_p = transpose(sub2ind(size(RGB),b_x_2,b_y_2));
end

function [RGB_marked] = markRGBwithFandBPoints(RGB, f_indices, b_indices)
    RGB(f_indices) = 256;
    RGB(b_indices) = 256;
    RGB_marked = RGB;
end

function [meets_requirements] = pixelMeetsRequirements(pixel_delta, Mask, r, c)
    % extract current pixel value    
    pixel_val = Mask(r,c);
    meets_requirements = true;
    
    % go through all boardering pixels with a radius of pixel_delta to
    % make sure its well centered
    for x = (r-pixel_delta):(r+pixel_delta)
        if x > 0 && x <= size(Mask,1)
            for y = (c-pixel_delta):(c+pixel_delta)
                if y > 0 && y <= size(Mask,2)
                    if Mask(x,y) ~= pixel_val
                        meets_requirements = false;
                    end
                end
            end
        end
    end
end

function [Mask] = calcMask(RGB_x, f_indices, b_indices)
    L = superpixels(RGB_x,500);
    Mask = lazysnapping(RGB_x,L,f_indices,b_indices);
end