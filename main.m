close all
clc
clear



%% Constants Settings
numPic = 18;
desired_numOfInterestPt = 1000;
nearestneighbor_threshold = 0.55;

% pic_start = 259;
% outputfile_str = sprintf('%d~%d', pic_start, numPic+pic_start-1);

%% Read in all images
image = cell(numPic,1);
for i=1:numPic
    filename = sprintf('image/prtn%02d.jpg', i-1);
    image{i} = imread(filename);
%     image{i} = imresize(image{i},0.5);
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
    filename = sprintf('image/proj/prtn%02d.jpg', i-1);
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

for img_num=numPic:-1:2
    
    descriptor1 = all_img_descriptors{img_num};
    descriptor2 = all_img_descriptors{img_num-1};
    
    nn_threshold = nearestneighbor_threshold;
    matches = featureMatching(descriptor1, descriptor2, nn_threshold);
    matrix = cell2mat(matches);
    img1_matched_points = matrix(1:end, 1:2);
    img2_matched_points = matrix(1:end, 3:4);
    
    % Save to cell for future processing
    matched_pairs{img_num, img_num-1} = {img1_matched_points;img2_matched_points};
    matched_pairs{img_num-1, img_num} = {img2_matched_points;img1_matched_points};
   
end

disp('Feature Matched Result');
toc
disp('---');

%% showMatchedFeatures 
% disp('Show Matched features...');
% tic
% 
%     for img_num=numPic:-1:2
%         match = matched_pairs{img_num, img_num-1};
%         
%         matched_points1 = match{1}(1:end, 2:-1:1);
%         matched_points2 = match{2}(1:end, 2:-1:1);
%     
%         figure;
%         showMatchedFeatures( image{img_num}, image{img_num-1}, matched_points1, matched_points2, 'montage');
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




%% Combine & Blending






%% Function: extractFeatures
function [feature_descriptor] = extractFeatures(image, interestpoint_num)
    % interestpoint_num: Desired number of interest points
    
    [fHMFirst, feature_record] = multiScaleHarris(image);
    
    [anms_record] = adaptiveNonMaximalSuppression(fHMFirst, feature_record, interestpoint_num);
    
    [feature_descriptor] = findMsopDescriptor(image, anms_record);
end
