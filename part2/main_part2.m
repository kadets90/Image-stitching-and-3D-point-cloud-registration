clear

camara_parameters = load('CalibrationData.mat');
%folder = 'both\board1';
%folder = 'both\board2';
%folder = 'both\room1';
%folder = 'both\table';
%folder = 'both\cenas';
folder = 'reconstruction\labpiv\data_rgb';
files = dir(folder);


cam_params.Kdepth = camara_parameters.Depth_cam.K;  % the 3x3 matrix for the intrinsic parameters for depth
cam_params.Krgb = camara_parameters.RGB_cam.K; % the 3x3 matrix for the intrinsic parameters for rgb
cam_params.R = camara_parameters.R_d_to_rgb; % the Rotation matrix from depth to RGB (extrinsic params)
cam_params.T = camara_parameters.T_d_to_rgb; % The translation from depth to RGB

number_images = size(files,1)/2 - 1;
files(1:2) = [];

imglistrgb = cell(number_images,1);
imglistdepth = cell(number_images,1);
for i = 1:number_images
     % get files names
    image_name = files(i+number_images).name;
    depth_name = files(i).name;
        
    % write the name of the files with directory
    imglistrgb{i} = [folder '\' image_name];
    imglistdepth{i} = [folder '\' depth_name];
end

[transforms, objects] = part2(imglistdepth, imglistrgb, cam_params);