%PIV first task - Matlab suggestion using SURF
clear

%% Load images.
%files_names = {'IMG_1786-2.jpg', 'IMG_1787-2.jpg', 'IMG_1788-2.jpg'};
%files_names = {'parede1.jpg', 'parede2.jpg'};
%files_names = {'sala1.jpg', 'sala2.jpg'};
%files_names = {'board1.jpg', 'board2.jpg'};
%files_names = {'pc1.jpg', 'pc2.jpg'};
%files_names = {'q1.jpg', 'q2.jpg', 'q3.jpg'};
%files_names = {'rua3.jpg', 'rua2.jpg', 'rua1.jpg'};
%files_names = {'BD1.PNG', 'BD2.PNG'};
%files_names = {'sintetico/im0001.jpg' 'sintetico/im0002.jpg' 'sintetico/im0003.jpg' 'sintetico/im0004.jpg' 'sintetico/im0005.jpg' 'sintetico/im0006.jpg' 'sintetico/im0007.jpg' 'sintetico/im0008.jpg' 'sintetico/im0009.jpg' 'sintetico/im0010.jpg'};
%files_names = {'sinteticoprojective/sinteticoprojective0001.jpg' 'sinteticoprojective/sinteticoprojective0002.jpg' 'sinteticoprojective/sinteticoprojective0003.jpg' 'sinteticoprojective/sinteticoprojective0004.jpg' 'sinteticoprojective/sinteticoprojective0005.jpg' 'sinteticoprojective/sinteticoprojective0006.jpg' 'sinteticoprojective/sinteticoprojective0007.jpg' 'sinteticoprojective/sinteticoprojective0008.jpg' 'sinteticoprojective/sinteticoprojective0009.jpg'};
%files_names = {'motion1af/IMG_20191007_133341.jpg'  'motion1af/IMG_20191007_133346.jpg'  'motion1af/IMG_20191007_133352.jpg'  'motion1af/IMG_20191007_133406.jpg'  'motion1af/IMG_20191007_133344.jpg'  'motion1af/IMG_20191007_133349.jpg'  'motion1af/IMG_20191007_133354.jpg'};
%files_names = {'motion2pro/IMG_20191007_133434.jpg'  'motion2pro/IMG_20191007_133437.jpg'  'motion2pro/IMG_20191007_133440.jpg'  'motion2pro/IMG_20191007_133442.jpg'  'motion2pro/IMG_20191007_133446.jpg'};
%files_names = {'translation/im_0001.jpg'  'translation/im_0002.jpg'  'translation/im_0003.jpg'  'translation/im_0004.jpg'  'translation/im_0005.jpg'  'translation/im_0006.jpg'  'translation/im_0007.jpg'  'translation/im_0008.jpg'  'translation/im_0009.jpg'  'translation/im_0010.jpg'};
%files_names = {'cmviana/IMG_20191019_145711.jpg' 'cmviana/IMG_20191019_145713.jpg' 'cmviana/IMG_20191019_145714.jpg' 'cmviana/IMG_20191019_145716.jpg' 'cmviana/IMG_20191019_145718.jpg' 'cmviana/IMG_20191019_145719.jpg' 'cmviana/IMG_20191019_145721.jpg' 'cmviana/IMG_20191019_145722.jpg' 'cmviana/IMG_20191019_145724.jpg'};
%files_names = {'vianaPiv\IMG_20191019_152859.jpg' 'vianaPiv\IMG_20191019_152858.jpg' 'vianaPiv\IMG_20191019_152857.jpg' 'vianaPiv\IMG_20191019_152856.jpg' 'vianaPiv\IMG_20191019_152855.jpg' 'vianaPiv\IMG_20191019_152854.jpg'};
%files_names = {'prepublica\IMG_20191019_145415.jpg' 'prepublica\IMG_20191019_145417.jpg' 'prepublica\IMG_20191019_145418.jpg' 'prepublica\IMG_20191019_145419.jpg' 'prepublica\IMG_20191019_145420.jpg' 'prepublica\IMG_20191019_145421.jpg' 'prepublica\IMG_20191019_145422.jpg' 'prepublica\IMG_20191019_145423.jpg' 'prepublica\IMG_20191019_145424.jpg' 'prepublica\IMG_20191019_145425.jpg'};

%% getting features
I = imread(char(files_names(1)));

% Initialize features for I(1)
grayImage1 = rgb2gray(I);
points = detectSURFFeatures(grayImage1);
[features, points] = extractFeatures(grayImage1, points);

% Initialize all the transforms to the identity matrix. Note that the
% projective transform is used here because the building images are fairly
% close to the camera. Had the scene been captured from a further distance,
% an affine transform would suffice.
numImages = size(files_names,2);

match_list = cell(numImages - 1, 2);

% Initialize variable to hold image sizes.
imageSize = zeros(numImages,2);

% Iterate over remaining image pairs
for n = 2:numImages
    
    % Store points and features for I(n-1).
    pointsPrevious = points;
    featuresPrevious = features;
        
    % Read I(n).
    I2 = imread(char(files_names(n)));
    
    % Convert image to grayscale.
    grayImage2 = rgb2gray(I2);    
    
    % Save image size.
    imageSize(n,:) = size(grayImage2);
    
    % Detect and extract SURF features for I(n).
    points = detectSURFFeatures(grayImage2);    
    [features, points] = extractFeatures(grayImage2, points);
  
    % Find correspondences between I(n) and I(n-1).
    indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);
       
    matchedPoints = points(indexPairs(:,1), :);
    matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);        
 
    % Store the matched points in a cell array
    match_list(n-1,1) = {matchedPointsPrev.Location};
    match_list(n-1,2) = {matchedPoints.Location};

end
%figure; ax = axes;
%showMatchedFeatures(grayImage1,grayImage2,matchedPointsPrev,matchedPoints,'montage','Parent',ax);

%% Our code
[H, world_pts_x, world_pts_y] = part1(files_names, match_list);

figure
plot([world_pts_x world_pts_x(:,1)]',[world_pts_y world_pts_y(:,1)]');
axis ij;axis equal;
% x = world_pts_x';
% x = x(:);
% y = world_pts_y';
% y = y(:);
% n_points = 4;
% plot(x(1:n_points), y(1:n_points), '*')
% hold on
% plot(x((n_points + 1):(n_points * 2)), y((n_points + 1):(n_points * 2)), '*')
% set(gca, 'YDir','reverse')
grid on
pause
%%
% Compute the output limits  for each transform

for i = 1:numImages
    tforms(i) = projective2d(H{i,1}');
end

for i = 1:numel(tforms)           
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);    
end


%%Make panaroma from center pixels
avgXLim = mean(xlim, 2);

[~, idx] = sort(avgXLim);

centerIdx = floor((numel(tforms)+1)/2);

centerImageIdx = idx(centerIdx);

Tinv = invert(tforms(centerImageIdx));

for i = 1:numel(tforms)    
    tforms(i).T = tforms(i).T * Tinv.T;
end

for i = 1:numel(tforms)           
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);
end

maxImageSize = max(imageSize);

% Find the minimum and maximum output limits 
xMin = min([1; xlim(:)]);
xMax = max([maxImageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([maxImageSize(1); ylim(:)]);

% Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);

% Initialize the "empty" panorama.
panorama = zeros([height width 3], 'like', I);

blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');  

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);

% Create the panorama.
for i = 1:numImages
    
    I = imread(char(files_names(i)));   
   
    % Transform I into the panorama.
    warpedImage = imwarp(I, tforms(i), 'OutputView', panoramaView);
                  
    % Generate a binary mask.    
    mask = imwarp(true(size(I,1),size(I,2)), tforms(i), 'OutputView', panoramaView);
    
    % Overlay the warpedImage onto the panorama.
    panorama = step(blender, panorama, warpedImage, mask);
end

figure
imshow(panorama)