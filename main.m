close all
% clc
% clear
% 
% %% Constants Settings
% numPic = 18;
% desired_numOfInterestPt = 1000;
% nearestneighbor_threshold = 0.7;
% 
% %% Read in all images
% image = cell(numPic,1);
% for i=1:numPic
%    
% %     resized_filename = sprintf('imageset/%d.jpg', i);
% %     if exist(resized_filename, 'file')
% %         image{i} = imread(resized_filename);
% %         continue;
% %     end
% 
%     filename = sprintf('imageset/%d.jpg', i); % image/prtn%02d.jpg
%     image{i} = imread(filename);
% %     image{i} = imresize(image{i},0.15);
% %     image{i} = rot90(image{i}, -1);
% %     imwrite(image{i}, resized_filename);
% %     imgInfo = imfinfo(filename)
% %     exposure_times(i) = imgInfo.DigitalCamera.ExposureTime;
% end
% % Test photo set
% % focal_length = [704.916; 706.286; 705.849; 706.645; 
% %                 706.587; 705.645; 705.327; 704.696;
% %                 703.794; 704.325; 704.696; 703.895;
% %                 704.289; 704.676; 704.847; 704.537;
% %                 705.102; 705.576];
% 
% focal_length = [
% 508.296 % ori11
% 506.089
% 502.303
% 505.51
% 506.93
% 
% 510.676
% 509.601
% 510.008
% 508.219
% 505.949
% % 512.234 %ori 1
% 506.66 
% 506.603
% 505.323
% 507.574 
% 
% 506.484
% 509.637
% 509.382
% 510.168
% % 507.529
% 
% ];
% %% Project to Cylindrical coordinate
% tic
% disp('Transform image to Cylindrical coordinate..');
% for i=1:numPic
%     filename = sprintf('imageset/proj/%d.jpg', i);
%     if exist(filename, 'file')
%         image{i} = imread(filename);
%     
%     else
%         image{i} = cylindricalProjection(image{i}, focal_length(i));
%         
%         imwrite(image{i}, filename);
%     end
%     
% end
% disp('All images Transformed');
% toc
% disp('---');
% 
% %% Extract all interest points & feature descriptors
% tic
% disp('Start Feature Extraction...');
% all_img_descriptors = cell(numPic,1);
% numberOfInterestPoint = desired_numOfInterestPt;
% 
% for img_num=1:numPic
%     all_img_descriptors{img_num} = extractFeatures(image{img_num}, numberOfInterestPoint);
%     fprintf('Processed image %d / %d ...\n', img_num, numPic);
% end
% disp('All features Extracted');
% toc
% disp('---');
% 
% %% Feature Matching & showMatchedFeatures
% 
% disp('Matching features...');
% tic
% 
% matched_pairs = cell(numPic, numPic);
% 
% for img_num=1:numPic-1
%     
%    
%     descriptor1 = all_img_descriptors{img_num};
%     descriptor2 = all_img_descriptors{img_num+1};
%     
%     nn_threshold = nearestneighbor_threshold;
%     matches = featureMatching(descriptor1, descriptor2, nn_threshold);
%     matrix = cell2mat(matches);
%     img1_matched_points = matrix(1:end, 1:2);
%     img2_matched_points = matrix(1:end, 3:4);
%     
%     % Save to cell for future processing
%     matched_pairs{img_num, img_num+1} = {img1_matched_points;img2_matched_points};
%     matched_pairs{img_num+1, img_num} = {img2_matched_points;img1_matched_points};
%    
% end
% 
% disp('Feature Matched Result');
% toc
% disp('---');
% 
% %% showMatchedFeatures 
% disp('Show Matched features...');
% tic
% 
%     for img_num=1:numPic-1
% 
%         match = matched_pairs{img_num, img_num+1};
%         
%         matched_points1 = match{1}(1:end, 2:-1:1);
%         matched_points2 = match{2}(1:end, 2:-1:1);
%     
%         figure;
%         showMatchedFeatures( image{img_num}, image{img_num+1}, matched_points1, matched_points2, 'montage');
%         
%     end
% disp('Matched Result Figures');
% toc
% disp('---');

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
    inlier_range = 80;
    for img_num=1:numPic-1
        
     
        match = matched_pairs{img_num, img_num+1}; 
        
        matched_points1 = match{1};
        matched_points2 = match{2};
        ransac_sample_times = round(size(matched_points1,1)/1.5);
        [mx,my] = findTranslationWithRansac(matched_points1, matched_points2, ransac_sample_times, inlier_range);
        translate(img_num, :) = [mx,my];
        
    end
    
disp('Finished Image Matching');
toc
disp('---');


%% Combine & Blending
previous = zeros(0,0,3);
translateVal_x_sum = 0;
translateVal_y_sum = 0;

for img_num=1:numPic-1
   
    img1 = image{img_num};
    if size(previous,1) ~= 0
        img1 = previous;
    end
    img2 = image{img_num+1};
    ROWS = size(img1,1);
    COLS = size(img1,2);
    ROWS_IMG2 = size(img2,1);
    COLS_IMG2 = size(img2,2);
    translateVal = translate(img_num,:);
    translateVal_x_sum = translateVal_x_sum+translateVal(1);
    translateVal_y_sum = translateVal_y_sum+translateVal(2);
   
    
    % extend x-dir(col)
    img1 = [img1, zeros(ROWS, translateVal(1) , 3)]; 
    
    img2 = [zeros(ROWS-ROWS_IMG2, COLS_IMG2 , 3); img2]; % Make row size same first
    img2 = [zeros(ROWS, translateVal_x_sum , 3), img2];
    
    % extend y-dir(row)
    if translateVal_y_sum >= 0
        img1 = [img1;zeros(translateVal_y_sum,size(img1,2),3)];
        img2 = [zeros(translateVal_y_sum,size(img2,2),3);img2];
    else
        img1 = [zeros(-translateVal_y_sum,size(img1,2),3);img1];
        img2 = [img2;zeros(-translateVal_y_sum,size(img2,2),3)];
    end
   
    extended_rowsize = size(img1,1);
    extended_colsize = size(img1,2);
    
    % Overlap Section: Ignore black frame
    window = 15;
    overlap1 = img1(:, translateVal_x_sum+window:COLS-window, :);

    % Calculate mask for linear blending
    overlap_width = size(overlap1,2);
    mask = zeros(extended_rowsize, overlap_width);
    mask2 = zeros(extended_rowsize, overlap_width);
    
    for i=1:overlap_width
        
        mask(:,i) = ones(extended_rowsize,1) .* (overlap_width-i)/(overlap_width-1.0);
        mask2(:,i) = ones(extended_rowsize,1) .* (i-1.0)/(overlap_width-1.0);
      
    end
    
    non_overlap_width1 = translateVal_x_sum+window;
    non_overlap_width2 = translateVal(1)+window-1;
    mask = [ones(extended_rowsize, non_overlap_width1), mask, zeros(extended_rowsize, non_overlap_width2)];
    mask2 = [zeros(extended_rowsize, non_overlap_width1), mask2, ones(extended_rowsize, non_overlap_width2)];

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


%% Crop Panoroma & Save
resultRowSize=  size(result,1);
resultColSize=  size(result,2);

result = result(resultRowSize-ROWS_IMG2+25:resultRowSize-abs(translateVal_y_sum)-25,15:resultColSize-15,:);

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
