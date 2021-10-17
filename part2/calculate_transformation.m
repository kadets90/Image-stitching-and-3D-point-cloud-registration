function T = calculate_transformation(X_1,Y_1,Z_1,X_2,Y_2,Z_2)
%with ransac, N random points are chosen, the centroid is computed in both
%frames, the transformation is computed
%we want the transformation from Y to X, X = i, Y = i+1, X = T*Y


n_matched_points = length(X_1);
N = 3;  % number of points randomly selected to calculate transformation
iterations = 100;
max_error = 0.01; % max distance from calculated point and real point
score = zeros(1,iterations);
score_threshold = 0.75; % if score > score_threshold, a good transformation was found
matching_points = 0;
best = -Inf;
index = 1;
T_vector = zeros(4,4,iterations); % vector of transformations for each iteration

for i=1:iterations
    
    k = randi(n_matched_points, 1, N);
    centroid1 = 1/N * [sum(X_1(k)); sum(Y_1(k)); sum(Z_1(k))];
    centroid2 = 1/N * [sum(X_2(k)); sum(Y_2(k)); sum(Z_2(k))];
    
    Xi = zeros(3,N);
    Yi = zeros(3,N);
    
    Xi(1,:) = X_1(k); Xi(2,:) = Y_1(k); Xi(3,:) = Z_1(k);
    Yi(1,:) = X_2(k); Yi(2,:) = Y_2(k); Yi(3,:) = Z_2(k);
    
    X = Xi - centroid1;
    Y = Yi - centroid2;
    
    S = Y*X';
    [U,E,V] = svd(S);
    Diag = eye(3,3);
    Diag(3,3) = det(V*U');
    R = V*Diag*U';
    
    translation = centroid1 - R*centroid2;
        
    T_return = zeros(4,4);
    T_return(1:3,1:3) = R(:,:)';
    T_return(4,1:3) = translation(:);
    T_return(4,4) = 1;
    
    % X = T*Y
    for k=1:length(X_1)
        dif = [X_1(k) Y_1(k) Z_1(k) 1] - [X_2(k) Y_2(k) Z_2(k) 1]*T_return;
        if norm(dif) < max_error
            % point matches
            matching_points = matching_points+1;
        end
    end
    
    score(i) = matching_points/length(X_1);
    
    if(score(i) >= score_threshold)
        T_vector(:,:,i) = T_return;
        index = i;
        break;
    elseif(score(i) > best)
        best = score(i);
        index = i;
        matching_points = 0;
    else
        matching_points = 0;
    end 

    T_vector(:,:,i) = T_return;
end

T = T_vector(:,:,index);
end

