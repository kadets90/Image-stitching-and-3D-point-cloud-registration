function [transforms, objects] = part2(imglistdepth, imglistrgb, cam_params)
    warning('off','all')
    number_images = max(size(imglistdepth));
    scale_factor = 1000; % jpg files use 1, for the png use 1000, converts mm to meters
    total_error = 0;
    total_inliers = 0;
    %% Find match_points
    rgb_image = imread(imglistrgb{1});
    depth_array = load(imglistdepth{1});
    depth_array = depth_array.depth_array;
    depth_array(isnan(depth_array)) = 0;

    % Initialize features for I(1)
    grayImage = rgb2gray(rgb_image);
    points = detectSURFFeatures(grayImage);
    [features, points] = extractFeatures(grayImage, points);

    % Calculating the rotation matrix - R - between every two point clouds
    % Calculating the Translation vector - T
    % Construct the transformation matrix - Transformation [R T; 0 0 0 1]
    last_line = [0 0 0 1];
    Transformation = cell(number_images-1,1);
    transforms = cell(number_images-1,1);
    R = cell(number_images-1,1);
    T = cell(number_images-1,1);

    % Iterate over remaining image pairs
    for n = 2:number_images

        % Store points and features for I(n-1).
        pointsPrevious = points;
        featuresPrevious = features;
        rgb_imagePrevious = rgb_image;
        depth_arrayPrevious = depth_array;

        % Read I(n).
        rgb_image = imread(imglistrgb{n});
        depth_array = load(imglistdepth{n});
        depth_array = depth_array.depth_array;
        depth_array(isnan(depth_array)) = 0;

        % Convert image to grayscale.
        grayImage = rgb2gray(rgb_image);    

        % Detect and extract SURF features for I(n).
        points = detectSURFFeatures(grayImage);    
        [features, points] = extractFeatures(grayImage, points);

        % Find correspondences between I(n) and I(n-1).
        indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);

        matchedPoints = points(indexPairs(:,1), :);
        matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);    

        %figure; ax = axes;
        %showMatchedFeatures(rgb_imagePrevious,rgb_image,matchedPointsPrev,matchedPoints,'montage','Parent',ax);
        
        %%%%%%%%%%%%%%%%%%%%%%
        % POSSIBLE MAT ERROR %
        %%%%%%%%%%%%%%%%%%%%%%
        idxPrevious = round(matchedPointsPrev.Location);
        idx = round(matchedPoints.Location);
        
        pcPrevious = get_3D_coordinates(idxPrevious, depth_arrayPrevious, cam_params, scale_factor);
        pc = get_3D_coordinates(idx, depth_array, cam_params, scale_factor);
        pt_ = pc;
        pt = pcPrevious;
        [pt, pt_, inliers] = ransac(pt, pt_);
        if pt == 0
            return
        end
        %disp(['Number of inliers from ', num2str(n), ' to ' , num2str(n-1), ': ' num2str(inliers)]);
        
        % Transformation between images
        i = n - 1;
        p = pt_;
        p_ = pt;

        com = ones(size(p,1),1)*mean(p,1);
        com_ = ones(size(p_,1),1)*mean(p_,1);

        % 3-by-n matrices
        % [x1 .... xn
        %  y1 .... yn
        %  z1 .... zn];
        A = (p_ - com_)';
        B = (p - com)';

        %SVD - (A*B' = U*sigma*V', R=UV')
        X = A*B';
        [U,~,V] = svd(X);

        r = U*V';
        % make sure that the determinant of R is equal to 1
        D = [1 0 0; 0 1 0; 0 0 det(r)];
        R{i} = U*D*V';
        T{i} = com_(1,:)' - R{i}*com(1,:)';

        Transformation{i} = [R{i} T{i}; last_line];
        transforms{i}.R = R{i};
        transforms{i}.T = T{i};

        
        error = [p_';ones(1,inliers)] - Transformation{i}*[p';ones(1,inliers)];
        total_error = total_error + norm(error);
        total_inliers = total_inliers + inliers;
        %disp(['Média do erro: ', num2str(norm(error)/inliers*1000)])
        
        
        if i > 1
            % T(1,n) = T(1,n-1) * T(n-1,n)
            Transformation{i} = Transformation{i-1}*Transformation{i};
            transforms{i}.R = Transformation{i}(1:3,1:3);
            transforms{i}.T = Transformation{i}(1:3,4);
        end
                
    end
    
    %disp(['Média do erro: ', num2str(total_error/total_inliers*1000)])
    
    %% Calculating the point clouds - pc - for each rgb/depth image
    pc = cell(number_images,1);
    for i = 1:number_images
        % reading files
        rgb_image = imread(imglistrgb{i});
        depth_array = load(imglistdepth{i});
        depth_array = depth_array.depth_array;
        depth_array(isnan(depth_array)) = 0;

        % getting cloud points
        pc{i} = get_point_cloud(rgb_image, depth_array, cam_params, scale_factor);
        %figure
        %showPointCloud(pc{i})
    end


    %% Transforming all points to same coordinate system

    transformed_points = cell(number_images-1,1);
    transformed_points{1} = pc{1}.Location';

    sp = size(pc{1}.Location); % size of the list of points
    xyzPoints = zeros(sp(1)*number_images, sp(2));
    colorPoints = zeros(sp(1)*number_images, sp(2), 'uint8');
    xyzPoints(1:sp(1),:) = transformed_points{1}';
    colorPoints(1:sp(1),:) = pc{1}.Color;
    for i = 2:number_images

        transformed_points{i} = Transformation{i-1}(1:3,1:3) * pc{i}.Location' + Transformation{i-1}(1:3,4)*ones(1,sp(1));
        xyzPoints((i-1)*sp(1) + 1 : i*sp(1),:) = transformed_points{i}';
        colorPoints((i-1)*sp(1) + 1 : i*sp(1),:) = pc{i}.Color;

    end

    ptCloud = pointCloud(xyzPoints,'Color',colorPoints);
    %figure
    %showPointCloud(ptCloud)

    %transforms = Transformation;
    objects = [];
end
