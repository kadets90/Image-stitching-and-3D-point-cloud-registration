%% Compute Homographies INVERTED
    
    % Least squares
    for i = 1:(number_images - 1)
        % matchpoints image 1
        x = match_points{i,2}(:,1);
        y = match_points{i,2}(:,2);
        
        % matchpoints image 2
        x_ = match_points{i,1}(:,1);
        y_ = match_points{i,1}(:,2);
        
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
            p_error = Y_{m} - X_{m}*beta{m};
            p_error = vec2mat(p_error,2);
            p_error = p_error(:,1).^2 + p_error(:,2).^2;
            error(m) = mean(p_error)*max(max_inliers)/max_inliers(m);
        end
        
        
        % compute homographies
        while min(error) < 10^5
            [~, indice] = min(error);
            switch indice
                case 1
                    H{i+1,i} = [beta{1}(1) beta{1}(2) beta{1}(3);
                                beta{1}(4) beta{1}(5) beta{1}(6);
                                beta{1}(7) beta{1}(8) 1        ];
                    if min(abs(eig(H{i+1,i}))) < 0.05
                        disp('Projective error');
                        error(indice) = 10^5;
                    else
                        disp('Projective');
                        break
                    end

                case 2
                    H{i+1,i} = [beta{2}(1) beta{2}(2) beta{2}(3);
                                beta{2}(4) beta{2}(5) beta{2}(6);
                                0          0          1        ];
                    if min(abs(eig(H{i+1,i}))) < 0.05
                        error(indice) = 10^5;
                        disp('Affine error');
                    else
                        disp('Affine');
                        break
                    end        

                case 3
                    H{i+1,i} = [beta{3}(1) -beta{3}(2) beta{3}(3);
                                beta{3}(2) beta{3}(1)  beta{3}(4);
                                0          0           1        ];
                    if min(abs(eig(H{i+1,i}))) < 0.05
                        error(indice) = 10^5;
                        disp('Euclidean similarity error');
                    else
                        disp('Euclidean similarity');
                        break
                    end
            end
        end
       
        
        for j = 1:i-1
            H{i+1,i-j} = H{i,i-j}*H{i+1,i};
        end

    end
    
    
    
    
    
    
                    A = [X_ -Y_];
                    [evec,~] = eig(A'*A);
                    H_ = reshape(evec(:,1),[3,3])';
                    H_ = H_/H_(end); % make H(3,3) = 1