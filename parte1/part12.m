function [H, world_pts_x, world_pts_y] = part12(image_list, match_list, varargin)    %K, pointsw


    percentage = 40;

    warning('off','all')

    if nargin ~= 2 && nargin ~= 4
        fprintf('Total number of inputs should be 2 or 4, and it is = %d\n', 2+length(varargin));
        H = 0;            
        world_pts_x = 0;  
        world_pts_y = 0; 
        return
    end
    
    number_images = size(image_list,2);
    
    %% Find match_points
    match_points = cell(number_images, number_images,2);
    max_points = 0;
    for m = 1:number_images
    IPrevious = imread(char(image_list{m}));

    % Initialize features for I(1)
    grayImage = rgb2gray(IPrevious);
    pointsPrevious = detectSURFFeatures(grayImage);
    [featuresPrevious, pointsPrevious] = extractFeatures(grayImage, pointsPrevious);

        % Iterate over remaining image pairs
        for n = m+1:number_images

            % Read I(n).
            I = imread(char(image_list{n}));

            % Convert image to grayscale.
            grayImage = rgb2gray(I);    

            % Detect and extract SURF features for I(n).
            points = detectSURFFeatures(grayImage);    
            [features, points] = extractFeatures(grayImage, points);

            % Find correspondences between I(n) and I(n-1).
            indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);

            matchedPoints = points(indexPairs(:,1), :);
            matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);        

            % Store the matched points in a cell array
            match_points(m,n,1) = {matchedPointsPrev.Location};
            match_points(m,n,2) = {matchedPoints.Location};
            max_points = max(max_points, max(size(matchedPoints.Location)));
            %figure; ax = axes;
            %showMatchedFeatures(IPrevious,I,matchedPointsPrev,matchedPoints,'montage','Parent',ax);
        end
    end
    
    %% Compute Homographies
    H = cell(number_images, number_images);
    % Least squares
    for i = 1:number_images - 1
        k = i+1;
        % matchpoints image 1
        x = match_points{i,k,1}(:,1);
        y = match_points{i,k,1}(:,2);

        % matchpoints image 2
        x_ = match_points{i,k,2}(:,1);
        y_ = match_points{i,k,2}(:,2);

        % getting columns
        zero = zeros(size(x));
        one = ones(size(x));

        if sum(one) < 4
            disp('Not enougth points')
            continue
        end

        % column 1
        c1 = [x zero]';
        c1 = c1(:);
        % column 2
        c2 = [y zero]';
        c2 = c2(:);
        % column 3
        c3 = [one zero]';
        c3 = c3(:);
        % column 4
        c4 = [zero x]';
        c4 = c4(:);
        % column 5
        c5 = [zero y]';
        c5 = c5(:);
        % column 6
        c6 = [zero one]';
        c6 = c6(:);
        % column 7
        c7 = [x.*x_ x.*y_]';
        c7 = -c7(:);
        % column 8
        c8 = [y.*x_ y.*y_]';
        c8 = -c8(:);

        X = [c1 c2 c3 c4 c5 c6 c7 c8];

        Y = [x_ y_]';
        Y = Y(:);

        % Applying RANSAC
        X_ = cell(3,1);
        Y_ = cell(3,1);
        n_pontos = [4 3 2];

        [X_{1}, Y_{1}, max_inliers(1)] = ransac(X, Y, n_pontos(1));
        [X_{2}, Y_{2}, max_inliers(2)] = ransac(X(:,1:6), Y, n_pontos(2));
        [X_{3}, Y_{3}, max_inliers(3)] = ransac([X(:,1) + X(:,5)  -X(:,2) + X(:,4)  X(:,3)  X(:,6)], Y, n_pontos(3));

        beta = cell(3,1);
        error = zeros(1,3);

        for m = 1:3
            beta{m} = (X_{m}'*X_{m})\(X_{m}'*Y_{m});

            %disp(['Inliers: ', num2str(max_inliers(m))])
            %eigen_values = eig(X_{m}'*X_{m});
            %disp(['Eigenvalues (max/min reason): ', num2str(sqrt(max(abs(eigen_values))/min(abs(eigen_values))))])
            %disp(['Minimum eigenvalue: ', num2str(min(abs(eigen_values)))])
            p_error = Y_{m} - X_{m}*beta{m};
            p_error = vec2mat(p_error,2);
            p_error = p_error(:,1).^2 + p_error(:,2).^2;
            %fprintf(['Errors\n - mean: ', num2str(mean(p_error)), '\n - min: ', num2str(min(p_error)), '\n - max: ', num2str(max(p_error)), '\n\n'])
            error(m) = mean(p_error)*max(max_inliers)/max_inliers(m);
        end


        % compute homographies
        H{i,i} = eye(3);
        while min(error) < 10^5
            [~, indice] = min(error);
            switch indice
                case 1
                    H{i,k} = [beta{1}(1) beta{1}(2) beta{1}(3);
                              beta{1}(4) beta{1}(5) beta{1}(6);
                              beta{1}(7) beta{1}(8) 1        ];
                    if min(abs(eig(H{i,k}))) < 0.05
                        disp('Projective error');
                        error(indice) = 10^5;
                    else
                        disp('Projective');
                        break
                    end

                case 2
                    H{i,k} = [beta{2}(1) beta{2}(2) beta{2}(3);
                                beta{2}(4) beta{2}(5) beta{2}(6);
                                0          0          1        ];
                    if min(abs(eig(H{i,k}))) < 0.05
                        error(indice) = 10^5;
                        disp('Affine error');
                    else
                        disp('Affine');
                        break
                    end        

                case 3
                    H{i,k} = [beta{3}(1) -beta{3}(2) beta{3}(3);
                                beta{3}(2) beta{3}(1)  beta{3}(4);
                                0          0           1        ];
                    if min(abs(eig(H{i,k}))) < 0.05
                        error(indice) = 10^5;
                        disp('Euclidean similarity error');
                    else
                        disp('Euclidean similarity');
                        break
                    end
            end
        end
        H{k,i} = inv(H{i,k});
    end
    H{i+1,i+1} = eye(3);

    
    %% Improve Homographies Results
    
    % Least squares
    for i = 1:number_images - 2
        for k = i+2:number_images
            % matchpoints image 1
            x = match_points{i,k,1}(:,1);
            y = match_points{i,k,1}(:,2);

            % matchpoints image 2
            x_ = match_points{i,k,2}(:,1);
            y_ = match_points{i,k,2}(:,2);
            
            
            % getting columns
            zero = zeros(size(x));
            one = ones(size(x));
            
            if sum(one) < percentage/100*max_points
                H{i,k} = H{k-1,k}*H{i,k-1}; 
                H{k,i} = inv(H{i,k});
                continue
            end
            
            if sum(one) < 4
                disp('Not enougth points')
                continue
            end

            % column 1
            c1 = [x zero]';
            c1 = c1(:);
            % column 2
            c2 = [y zero]';
            c2 = c2(:);
            % column 3
            c3 = [one zero]';
            c3 = c3(:);
            % column 4
            c4 = [zero x]';
            c4 = c4(:);
            % column 5
            c5 = [zero y]';
            c5 = c5(:);
            % column 6
            c6 = [zero one]';
            c6 = c6(:);
            % column 7
            c7 = [x.*x_ x.*y_]';
            c7 = -c7(:);
            % column 8
            c8 = [y.*x_ y.*y_]';
            c8 = -c8(:);

            X = [c1 c2 c3 c4 c5 c6 c7 c8];

            Y = [x_ y_]';
            Y = Y(:);

            % Applying RANSAC
            X_ = cell(3,1);
            Y_ = cell(3,1);
            n_pontos = [4 3 2];

            [X_{1}, Y_{1}, max_inliers(1)] = ransac(X, Y, n_pontos(1));
            [X_{2}, Y_{2}, max_inliers(2)] = ransac(X(:,1:6), Y, n_pontos(2));
            [X_{3}, Y_{3}, max_inliers(3)] = ransac([X(:,1) + X(:,5)  -X(:,2) + X(:,4)  X(:,3)  X(:,6)], Y, n_pontos(3));

            beta = cell(3,1);
            error = zeros(1,3);

            for m = 1:3
                beta{m} = (X_{m}'*X_{m})\(X_{m}'*Y_{m});

                %disp(['Inliers: ', num2str(max_inliers(m))])
                %eigen_values = eig(X_{m}'*X_{m});
                %disp(['Eigenvalues (max/min reason): ', num2str(sqrt(max(abs(eigen_values))/min(abs(eigen_values))))])
                %disp(['Minimum eigenvalue: ', num2str(min(abs(eigen_values)))])
                p_error = Y_{m} - X_{m}*beta{m};
                p_error = vec2mat(p_error,2);
                p_error = p_error(:,1).^2 + p_error(:,2).^2;
                %fprintf(['Errors\n - mean: ', num2str(mean(p_error)), '\n - min: ', num2str(min(p_error)), '\n - max: ', num2str(max(p_error)), '\n\n'])
                error(m) = mean(p_error)*max(max_inliers)/max_inliers(m);
            end


            % compute homographies
            H{i,i} = eye(3);
            while min(error) < 10^5
                [~, indice] = min(error);
                switch indice
                    case 1
                        H{i,k} = [beta{1}(1) beta{1}(2) beta{1}(3);
                                  beta{1}(4) beta{1}(5) beta{1}(6);
                                  beta{1}(7) beta{1}(8) 1        ];
                        if min(abs(eig(H{i,k}))) < 0.05
                            disp('Projective error');
                            error(indice) = 10^5;
                        else
                            disp('Projective');
                            break
                        end

                    case 2
                        H{i,k} = [beta{2}(1) beta{2}(2) beta{2}(3);
                                    beta{2}(4) beta{2}(5) beta{2}(6);
                                    0          0          1        ];
                        if min(abs(eig(H{i,k}))) < 0.05
                            error(indice) = 10^5;
                            disp('Affine error');
                        else
                            disp('Affine');
                            break
                        end        

                    case 3
                        H{i,k} = [beta{3}(1) -beta{3}(2) beta{3}(3);
                                    beta{3}(2) beta{3}(1)  beta{3}(4);
                                    0          0           1        ];
                        if min(abs(eig(H{i,k}))) < 0.05
                            error(indice) = 10^5;
                            disp('Euclidean similarity error');
                        else
                            disp('Euclidean similarity');
                            break
                        end
                end
            end
            H{k,i} = inv(H{i,k});
        end
    end
    
    %% CALCULATE WORLD COORDINATES 
    if nargin == 4
        % CALCULATE WORLD TRANSFORM 
        K = varargin{1};
        pointsw = varargin{2}';
        x_ = pointsw(1,:)';
        y_ = pointsw(2,:)';
        Yw = pointsw(:);
        
        I = imread(char(image_list(1)));   % Read I(n).
        grayImage = rgb2gray(I);           % Convert image to grayscale.
        coordinates = size(grayImage);     % Save image size.
        corners = [1 1;
                   coordinates(2) 1;
                   1 coordinates(1);
                   coordinates(2) coordinates(1)];
        
        x = corners(:,1);
        y = corners(:,2);
        
        % getting columns
        zero = zeros(size(x));
        one = ones(size(x));
        
        % column 1
        c1 = [x zero]';
        c1 = c1(:);
        % column 2
        c2 = [y zero]';
        c2 = c2(:);
        % column 3
        c3 = [one zero]';
        c3 = c3(:);
        % column 4
        c4 = [zero x]';
        c4 = c4(:);
        % column 5
        c5 = [zero y]';
        c5 = c5(:);
        % column 6
        c6 = [zero one]';
        c6 = c6(:);
        % column 7
        c7 = [x.*x_ x.*y_]';
        c7 = -c7(:);
        % column 8
        c8 = [y.*x_ y.*y_]';
        c8 = -c8(:);
        
        Xw = [c1 c2 c3 c4 c5 c6 c7 c8];
            
        betaw = (Xw'*Xw)\(Xw'*Yw);
        
        Hw = [betaw(1) betaw(2) betaw(3);
              betaw(4) betaw(5) betaw(6);
              betaw(7) betaw(8) 1      ];
    else
        K = eye(3);
        Hw = eye(3);
    end
    
    world_pts_x = zeros(number_images,4);
    world_pts_y = zeros(number_images,4);
    
    % computing coordinates
    for j = 1:number_images 
        I = imread(char(image_list(j)));   % Read I(n).
        grayImage = rgb2gray(I);           % Convert image to grayscale.
        coordinates = size(grayImage);     % Save image size.

        UL = Hw * H{j,1} * [1 1 1]';
        UR = Hw * H{j,1} * [coordinates(2) 1 1]';
        DL = Hw * H{j,1} * [1 coordinates(1) 1]';
        DR = Hw * H{j,1} * [coordinates(2) coordinates(1) 1]';
        corners = [DL DR UR UL];
        world_pts_x(j,:) = corners(1,:)./corners(3,:);
        world_pts_y(j,:) = corners(2,:)./corners(3,:);
    end

end