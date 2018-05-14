function [feature_window_cell] = findMsopDescriptor(image, feature_record)
    %% Find Feature descriptor
    % MSOP: Sample 8x8 oriented patch from 40x40 window
    image = rgb2gray(image);
    img_size = size(image);
    
    ROWS = img_size(1);
    COLS = img_size(2);
    smoothedImage = imgaussfilt( image, 4.5 );
    [Gmag,Gdir] = imgradient(smoothedImage);

    count = 0;
    for i=1:ROWS
        for j=1:COLS
            
            if feature_record(i,j) > 0
                squareWindowAroundFeature = zeros(40, 40);
                sampleFeatureWindow = zeros(8, 8);
                
                if_eliminate_feature = 0;
                
                
                dir_angle = (Gdir(i,j)/180.0) * pi;
                dx = cos(dir_angle); % unit length
                dy = sin(dir_angle);
                start_x = floor(j - 20*dx);
                start_y = floor(i - 20*dy);
                
                for m=1:40
                    for n=1:40
                        x = floor(start_x + dx*m);
                        y = floor(start_y + dy*m);
                        
                        if x<1 || x>COLS || y<1 || y>ROWS
                            squareWindowAroundFeature(n,m) = 0;
                            if_eliminate_feature = 1;
                            break;
                        end
%                         disp(y);disp(x);
                        squareWindowAroundFeature(n,m) = image(y,x);
                    end
                    if(if_eliminate_feature==1); break; end
                end
                if(if_eliminate_feature==1); break; end
                % Sample to 8x8 patch
                for m=1:8
                    for n=1:8
                        start_h = (m-1)*5 +1;
                        start_w = (n-1)*5 +1;
                        sub_window = squareWindowAroundFeature(start_h:start_h+4, start_w:start_w+4);
                        meanVal = (sum(sub_window(:)))/25.0;

                        sampleFeatureWindow(m,n) = meanVal;
                    end
                end
                count = count+1;
                feature_window_cell{count,1} = i;
                feature_window_cell{count,2} = j;
                feature_window_cell{count,3} = sampleFeatureWindow;
            end
        end
    end

end