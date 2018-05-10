close all
clc
clear



%% Constants Settings
numPic = 18;
desired_numOfInterestPt = 1000;
nearestneighbor_threshold = 0.6;

% pic_start = 259;
% outputfile_str = sprintf('%d~%d', pic_start, numPic+pic_start-1);

%% Read in all images
image = cell(numPic,1);
for i=1:numPic
    filename = sprintf('image/prtnn%02d.jpg', i); % image/prtn%02d.jpg
    image{i} = imread(filename);
%     image{i} = imresize(image{i},0.1);
%     image{i} = rot90(image{i}, -1);
%     imgInfo = imfinfo(filename);
%     exposure_times(i) = imgInfo.DigitalCamera.ExposureTime;
end
focal_length = [704.916; 706.286; 705.849; 706.645; 
                706.587; 705.645; 705.327; 704.696;
                703.794; 704.325; 704.696; 703.895;
                704.289; 704.676; 704.847; 704.537;
                705.102; 705.576];

%% Project to Cylindrical coordinate
tic
disp('Transform image to Cylindrical coordinate..');
for i=1:numPic
    filename = sprintf('image/proj/prtnn%02d.jpg', i);
    if exist(filename, 'file')
        image{i} = imread(filename);
    
    else
        image{i} = cylindricalProjection(image{i}, focal_length(i));
        imwrite(image{i}, filename);
    end
    
end
disp('All images Transformed');
toc
disp('---');

%% Extract all interest points & feature descriptors
tic
disp('Start Feature Extraction...');
all_img_descriptors = cell(numPic,1);
numberOfInterestPoint = desired_numOfInterestPt;

for img_num=1:numPic
    all_img_descriptors{img_num} = extractFeatures(image{img_num}, numberOfInterestPoint);
    fprintf('Processed image %d / %d ...\n', img_num, numPic);
end
disp('All features Extracted');
toc
disp('---');

%% Feature Matching & showMatchedFeatures

disp('Matching features...');
tic

matched_pairs = cell(numPic, numPic);

for img_num=1:numPic-1
    
    descriptor1 = all_img_descriptors{img_num};
    descriptor2 = all_img_descriptors{img_num+1};
    
    nn_threshold = nearestneighbor_threshold;
    matches = featureMatching(descriptor1, descriptor2, nn_threshold);
    matrix = cell2mat(matches);
    img1_matched_points = matrix(1:end, 1:2);
    img2_matched_points = matrix(1:end, 3:4);
    
    % Save to cell for future processing
    matched_pairs{img_num, img_num+1} = {img1_matched_points;img2_matched_points};
    matched_pairs{img_num+1, img_num} = {img2_matched_points;img1_matched_points};
   
end

disp('Feature Matched Result');
toc
disp('---');

%% showMatchedFeatures 
disp('Show Matched features...');
tic

    for img_num=1:numPic-1
        match = matched_pairs{img_num, img_num+1};
        
        matched_points1 = match{1}(1:end, 2:-1:1);
        matched_points2 = match{2}(1:end, 2:-1:1);
    
        figure;
        showMatchedFeatures( image{img_num}, image{img_num+1}, matched_points1, matched_points2, 'montage');
    end
disp('Matched Result Figures');
toc
disp('---');

%% Image Matching with RANSAC
    % Matched results stored in matched_pairs
    % With matched feature location(row,col) in both image
    % Usage:
    %   match = matched_pairs{img_num_1, img_num_2}; 
    %   matched_points1 = match{1}; --> nx2 matrix (row,col)
    %   matched_points2 = match{2}; --> nx2 matrix (row,col)
disp('Image Matching...');
tic
    translate = zeros(numPic,2);
    inlier_range = 100;
    for img_num=1:numPic-1
        match = matched_pairs{img_num, img_num+1}; 
        
        matched_points1 = match{1};
        matched_points2 = match{2};
        ransac_sample_times = round(size(matched_points1,1));
        [mx,my] = findTranslationWithRansac(matched_points1, matched_points2, ransac_sample_times, inlier_range);
        translate(img_num, :) = [mx,my];
    end
    
disp('Finished Image Matching');
toc
disp('---');


%% Combine & Blending
previous = zeros(0,0,3);
translateVal_sum = 0;
for m=1:numPic-1
    img1 = image{m};
    if size(previous,1) ~= 0
        img1 = previous;
    end
    img2 = image{m+1};
    ROWS = size(img1,1);
    COLS = size(img1,2);
    ROWS_IMG2 = size(img2,1);
    COLS_IMG2 = size(img2,2);
    translateVal = translate(m,:);
    translateVal_sum = translateVal_sum+translateVal(1);
    
    % extend x-dir(col)
    img1 = [img1, zeros(ROWS, translateVal(1) , 3)]; 
    
    img2 = [zeros(ROWS-ROWS_IMG2, COLS_IMG2 , 3); img2]; % Make row size same first
    img2 = [zeros(ROWS, translateVal_sum , 3), img2];
    
    % extend y-dir(row)
    if translateVal(2) >= 0
        img1 = [img1;zeros(translateVal(2),size(img1,2),3)];
        img2 = [zeros(translateVal(2),size(img2,2),3);img2];
    else
        img1 = [zeros(-translateVal(2),size(img1,2),3);img1];
        img2 = [img2;zeros(-translateVal(2),size(img2,2),3)];
    end
    
    extended_rowsize = size(img1,1);
    extended_colsize = size(img1,2);
    overlap1 = img1(:, translateVal_sum:COLS, :);
%     overlap2 = img2(:, translateVal(1):COLS, :);

    % Calculate mask for linear blending
    overlap_width = size(overlap1,2);
    imshow(overlap1);
    mask = zeros(extended_rowsize, overlap_width);
    mask2 = zeros(extended_rowsize, overlap_width);
    for i=1:overlap_width

        mask(:,i) = ones(extended_rowsize,1) .* (overlap_width-i)/(overlap_width-1.0);
        mask2(:,i) = ones(extended_rowsize,1) .* (i-1.0)/(overlap_width-1.0);
    end
    
    mask = [ones(extended_rowsize, translateVal_sum-1), mask, zeros(extended_rowsize, translateVal_sum)];
    mask2 = [zeros(extended_rowsize, translateVal_sum-1), mask2, ones(extended_rowsize, translateVal_sum)];


    % Combine two image
    result = zeros(extended_rowsize,extended_colsize, 3, 'uint8');
    for i=1:extended_rowsize
        for j=1:extended_colsize
            img_mask_pixel1(1) = mask(i,j) * img1(i,j,1);
            img_mask_pixel1(2) = mask(i,j) * img1(i,j,2);
            img_mask_pixel1(3) = mask(i,j) * img1(i,j,3);

            img_mask_pixel2(1) = mask2(i,j) * img2(i,j,1);
            img_mask_pixel2(2) = mask2(i,j) * img2(i,j,2);
            img_mask_pixel2(3) = mask2(i,j) * img2(i,j,3);

            result(i,j,1) = uint8(img_mask_pixel1(1) + img_mask_pixel2(1));
            result(i,j,2) = uint8(img_mask_pixel1(2) + img_mask_pixel2(2));
            result(i,j,3) = uint8(img_mask_pixel1(3) + img_mask_pixel2(3));

        end
    end
    
    previous = result;
    
end
figure;
imshow(result);
title('Full panoroma');
imwrite(result, 'result/panoroma.jpg');

%% Function: extractFeatures
function [feature_descriptor] = extractFeatures(image, interestpoint_num)
    % interestpoint_num: Desired number of interest points
    
    [fHMFirst, feature_record] = multiScaleHarris(image);
    
    [anms_record] = adaptiveNonMaximalSuppression(fHMFirst, feature_record, interestpoint_num);
    
    [feature_descriptor] = findMsopDescriptor(image, anms_record);
end
