function [first_fHM, feature_record] = multiScaleHarris(image)

    img_size = size(image);
    ROWS = img_size(1);
    COLS = img_size(2);
    image_gray = rgb2gray(image);

    % Use Sobel operator to do derivative
    sobel_y = [-1 -2 -1; 0 0 0; 1 2 1];
    sobel_x = transpose(sobel_y);

    grad_x = [-1 0 1];
    grad_y = transpose(grad_x);
    sub_image = image_gray;
    feature_record = zeros(ROWS,COLS);
    
    first_fHM = zeros(ROWS,COLS);
    
    %% Do multi scale Harris Corner
    for lv=1:4
        %% Sub image Size & scale
        scale = 2^(lv-1);
        sub_img_size = size(sub_image);

        SUB_ROWS = sub_img_size(1);
        SUB_COLS = sub_img_size(2);
        
        %% Gradient of x&y
        Ix = filter2(sobel_x, sub_image);
        Iy = filter2(sobel_y, sub_image);

        %% outer product of gradient for M=[A C; C B]
        Ix2 = Ix .* Ix;
        Iy2 = Iy .* Iy;
        Ixy = Ix .* Iy;

        %% apply gaussian filter
        sigma = 1.5;
        Sx2 = imgaussfilt(Ix2, sigma); 
        Sy2 = imgaussfilt(Iy2, sigma);
        Sxy = imgaussfilt(Ixy, sigma);

        %% Compute corner strength: Fhm
        for i=1:SUB_ROWS
            for j=1:SUB_COLS
                M = [Sx2(i,j) Sxy(i,j); Sxy(i,j) Sy2(i,j)];
                fHM(i,j) = det(M) ./ trace(M);
            end
        end
        
        if lv==1
            first_fHM = fHM;
        end

        %% Pick local maxima of 3x3 and above threshold t
        window = 1;
        t = 10.0; % threshold
        
        for i=1:SUB_ROWS
            for j=1:SUB_COLS
                % Check each pixel start
                if fHM(i,j) > t
                    picked = 1;
                    if i-window<1 || i+window > SUB_ROWS || j-window<1 || j+window > SUB_COLS
                        continue;
                    end

                    neighbour3x3 = fHM(i-window:i+window,j-window:j+window);
                    if max(neighbour3x3(:)) ~= fHM(i,j)
                        picked = 0;
                    end
           
                    if picked == 1
                        r = i*scale;
                        c = j*scale;
                        feature_record(r,c) = 1;
                    end

                end

                % Check each pixel END
            end
        end
        
        %% image pyramid
        sub_image = impyramid(sub_image, 'reduce');
    end

end