I1 = rgb2gray(imread('sintetico/im0001.jpg'));
I2 = rgb2gray(imread('sintetico/im0002.jpg'));

% Detect SURF features

points1 = detectSURFFeatures(I1);
points2 = detectSURFFeatures(I2);

% Extract features

[f1, vpts1] = extractFeatures(I1, points1);
[f2, vpts2] = extractFeatures(I2, points2);

% Match features.

indexPairs = matchFeatures(f1, f2, 'MaxRatio', 0.5);
matchedPoints1 = vpts1(indexPairs(:, 1));
matchedPoints2 = vpts2(indexPairs(:, 2));

% Visualize candidate matches.

figure; ax = axes;
showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2,'montage','Parent',ax);
title(ax, 'Candidate point matches');
legend(ax, 'Matched points 1','Matched points 2');