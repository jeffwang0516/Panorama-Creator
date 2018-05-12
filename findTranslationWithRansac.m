function [mx,my] = findTranslationWithRansac(matched_points1, matched_points2, sample_times_K, inlier_range)

    % Consider translation only
    % assume stitching from left to right
    matchPairSize = size(matched_points1,1);
    
    % Translation offset
    mx = -1;
    my = -1;
    
    final_inlier_count = 0;
    for iter=1:sample_times_K
        
        mx_iter = -1;
        my_iter = -100;
        while mx_iter < 0 
            random_index = floor(rand * matchPairSize) + 1;
            point1 = matched_points1(random_index,:);
            point2 = matched_points2(random_index,:);
            mx_iter = point1(2) - point2(2);
            my_iter = point1(1) - point2(1);
        end
        
        inlier_count = 0;
        for i=1:matchPairSize
            if i == random_index; continue; end
            
            % Test Moving the left image
            point1 = matched_points1(i,:);
            point2 = matched_points2(i,:);
            x = (point1(2) - mx_iter) - point2(2);
            y = (point1(1) - my_iter) - point2(1);
            
            diff = x^2 + y^2;
            
            if diff < inlier_range
                inlier_count = inlier_count + 1;
            end
  
        end
        
        if inlier_count > final_inlier_count
            mx = mx_iter;
            my = my_iter;
            final_inlier_count = inlier_count;
         
        end
    end
%     fprintf("mx=%d, my=%d\n", mx,my);
end