function [X__, Y__, max_inliers] = ransac(X, Y, n_points)
                
    % Applying RANSAC
    max_inliers = 0;
    threshold = 9; % 3 pixeis
    max_iterations = 100;
    
    
    number_of_points = size(X,1)/2;
    if number_of_points < n_points
            X__ = X;
            Y__ = Y;
            return
    end
    
    y = vec2mat(Y,2);
    
    for k = 1:max_iterations
        % generating a set of points
        r = randi([1 number_of_points],1,n_points);
        while size(unique(X(2*r-1,1:2),'rows'),1) ~= n_points || size(unique(y(r,:),'rows'),1) ~= n_points
            r = randi([1 number_of_points],1,n_points);
        end

        r_ = [2*r-1 2*r];
        
        X_ = X(r_,:);
        Y_ = Y(r_);
        
        beta_ = (X_'*X_)\(X_'*Y_);

        % removing outliers
        error = vec2mat(Y - X*beta_,2);
        error = error(:,1).^2 + error(:,2).^2;
        inliers = length(find(error < threshold));

        % save the best set
        if inliers > max_inliers
            max_inliers = inliers;
            error = [error error]';
            error = error(:);
            X__ = X(error < threshold,:);
            Y__ = Y(error < threshold);
        end

    end
        
end