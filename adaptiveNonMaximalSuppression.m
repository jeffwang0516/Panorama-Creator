function [feature_anms_selected] = adaptiveNonMaximalSuppression(fHM, feature_record, num_of_interestpoints)
    % Input
    % fHM: corner strength matrix
    % feature_record: mask of selected points (2D), 0 or 1 for each position (row,col)
    % num_of_interestpoints: number of interest points to find
    %
    % Output
    % feature_anms_selected: 2D, same as feature_record

    img_size = size(fHM);
    ROWS = img_size(1);
    COLS = img_size(2);
    
    %% Pick out fHM values of selected features, reshape to nx1
    feature_index = find(feature_record > 0);
    featureVal = fHM( feature_index );

    
    %% sort in descending order
    [sorted_featureVal, sorted_order] = sort(featureVal, 'descend');
    
    sorted_feature_index = feature_index(sorted_order);
    num_of_feature = size(featureVal,1);
    
    %% transform to 2D axis, x & y
    x = mod(sorted_feature_index, ROWS);
    x(x==0) = ROWS;
    y = ceil(sorted_feature_index / ROWS);

    %% find min suppress radius R(i) for each interest point i

    for i=num_of_feature:-1:1
        min_dist = inf;

        if i>1
            for j=i-1:-1:1
                if fHM(sorted_feature_index(i)) < 0.9*fHM(sorted_feature_index(j))

                    dist = sqrt((x(j)-x(i))^2 + (y(j)-y(i))^2);
                    if min_dist > dist
                        min_dist = dist;
                    end                   
                end
            end
        end
        R(i) = min_dist;
    end

    %% Select n=num_of_interestpoints points with largest R(i)
    radius = 0;
    
    numOfIp = num_of_feature;
    
    while numOfIp >= num_of_interestpoints 
        index_of_features = find(R >= radius);
        numOfIp = size(index_of_features,2);
        radius = radius+0.1;
    end
    
    %% Construct output mask for selected features
    feature_anms_index = sorted_feature_index(index_of_features);
    feature_anms_selected = zeros(ROWS,COLS);
    for i=1:numOfIp
        x1= (mod(feature_anms_index(i), ROWS));
        if x1==0
            x1=ROWS;
        end
        x2= (ceil(feature_anms_index(i) / ROWS));
        feature_anms_selected(x1,x2)=1;
    end
    
end