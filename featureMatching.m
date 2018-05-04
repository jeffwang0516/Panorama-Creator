function [match_pair] = featureMatching(descriptor1, descriptor2, nn_threshold)
    % Exhaustively calculate distance between a feature in one image to the other.
    % Thresholding on the ratio e1-nn & e2-nn
    % [Brown, Szeliski, Winder, CVPR'2005] 5.1 section
    
    descriptor1_size = size(descriptor1,1);
    descriptor2_size = size(descriptor2,1);
    
    count =1;
    dist = zeros(descriptor1_size,descriptor2_size);
    for i=1:descriptor1_size

        for j=1:descriptor2_size
            difference_mat = descriptor1{i,3} - descriptor2{j,3};
            
            dist(i,j) = norm(difference_mat(:));
        end        
    end
    
    %% nearest neighbor ratio check
    threshold = nn_threshold;
    
    diff = zeros(descriptor1_size,descriptor2_size);
    index_of_desc2 = zeros(descriptor1_size,descriptor2_size);
    geometric_dist = zeros(descriptor1_size);
    nearest_neighbor_ratio = zeros(descriptor1_size);
    
    for i=1:descriptor1_size
        
        [diff(i,:), index_of_desc2(i,:)] = sort(dist(i,:), 'ascend');
        [x1,y1] = descriptor1{i,1:2};
        [x2,y2] = descriptor2{index_of_desc2(i,1), 1:2};
%         geometric_dist(i) = norm([x1-x2, y1-y2]);
        geometric_dist(i) = sqrt((x1-x2)^2 + (y1-y2)^2);
        
        nearest_neighbor_ratio(i) = diff(i,1) / diff(i,2);
    end
    
    for i=1:descriptor1_size
        if nearest_neighbor_ratio(i) < threshold
            if diff(i,1) < 100
                match_pair{count,1} = [descriptor1{i,1:2}];
                match_pair{count,2} = [descriptor2{index_of_desc2(i,1),1:2}];
                match_pair{count,3} = diff(i,1);
                
                count = count+1;
            end
        end
    end
    
end

