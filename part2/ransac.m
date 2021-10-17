function [pt, pt_, max_inliers] = ransac(pc, pc_)


    % Applying RANSAC
    max_inliers = 0;
    threshold = 0.1; % 10 cm
    max_iterations = 100;
    n_points = 4;
    
    number_of_points = size(pc,1);
    if number_of_points < n_points
        disp('Not enough match points!')
        pt=0;
        pt_=0;
        return;
    end
    for k = 1:max_iterations
        % generating a set of points
        r = randi([1 number_of_points],1,n_points);

        while size(unique(pc(r,:),'rows'),1) ~= n_points || size(unique(pc_(r,:),'rows'),1) ~= n_points
            r = randi([1 number_of_points],1,n_points);
        end
        
        p = pc(r,:);
        p_ = pc_(r,:);
        
        % center of mass
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
        R = U*D*V';
        

        T = com_(1,:)' - R*com(1,:)';

        % removing outliers
        error = pc_' - (R*pc' + T*ones(1,number_of_points));
        error = error(1,:).^2 + error(2,:).^2 + error(3,:).^2;
        inliers = length(find(error < threshold^2));
%% CHECK NUMBER OF INLIERS IF < 3 SOMETHING IS WRONG
% PAY ATTENTION TO IMAGINARY VALUES

        % save the best set
        if inliers > max_inliers
            max_inliers = inliers;
            pt = pc((error < threshold^2),:);
            pt_ = pc_((error < threshold^2),:);
        end

    end
        
end