function [projs] = crop_images_box(path_box_file, path_micrograph, diameter)
%% KLTpicker_to_images
% Generate particle images from micrograph's picker box file.
% Function uses ASPIRE functions such as readSTAR, ReadMRC, etc.
% 
%   Input:
%            path_star_file:  full path to KLT-picker box file (*.box)
%            path_micrograph: full path to MRC micrograph (*.mrc) or
%                             an LxLx1 double variable
%            diameter:        of each cropped particle image (r_max * 2)
%
%   Output:  projs:           size LxLxN
% 
% Versions:
% 0.1        |  Gili Weiss-Dicker, April 2021 
%% Configurations
plot_micrograph        = 0;
plot_projs             = 1;

%% Read KLT picker results
if ~exist('path_box_file','var') || ~exist('path_micrograph','var')
    error('crop_images_box can''t run without box file or micrograph')
end
if ~exist('diameter','var')
    diameter=360; % 10028 data, must be an even number
end
log_message('Running crop_images_box...')

if ischar(path_micrograph)    % path to *.mrc input file 
    micrograph = ReadMRC(path_micrograph);
else                          % variable input
    micrograph = path_micrograph;
end

datablocks     = load(path_box_file);

rlnCoordinateX   = datablocks(:,1);
rlnCoordinateY   = datablocks(:,2);

invalid_x_inds = find (rlnCoordinateX< 1); 
invalid_y_inds = find (rlnCoordinateY< 1); 
invalid_oob_x_inds = find (rlnCoordinateX> (size(micrograph,1) - diameter)); 
invalid_oob_y_inds = find (rlnCoordinateY> (size(micrograph,1) - diameter));
all_inds = [ invalid_x_inds; invalid_y_inds;invalid_oob_x_inds; invalid_oob_y_inds];
rlnCoordinateX(all_inds)= [];
rlnCoordinateY(all_inds)= [];

if plot_micrograph
    figure; imshow(micrograph, []);
end

%% Crop particles from micrograph
Nprojs            = size(rlnCoordinateX, 1);
r                 = diameter/2;               % image radious
projs             = zeros( diameter, diameter, size(rlnCoordinateX, 2));

micrograph_fixed  = micrograph; %imrotate(micrograph,90);  % rotate (to match EMAN)
count = 0;
for ii = 1:Nprojs
    % particle position on micrograph
    center_position = [rlnCoordinateX(ii), rlnCoordinateY(ii)]; 

    % crop particle images
    projs(:,:,ii) = micrograph_fixed( center_position(1) : center_position(1)+2*r-1 , center_position(2): center_position(2)+2*r-1);
    count = count +1;
end

if plot_projs
    figure; viewstack(projs,10,10,0); title('Extracted particle images')
end
log_message('Finished image cropping from micrograph');
end