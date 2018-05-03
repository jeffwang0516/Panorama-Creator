function image_cylCoord = cylindricalProjection(image, focal_length)
    img_size = size(image);
    ROWS = img_size(1);
    COLS = img_size(2);
    
    image_cylCoord = zeros(img_size, 'uint8');
    
    s = double(focal_length);
  
    for i=1:ROWS
        for j=1:COLS
            % Let the center of image be (0,0)
            y = round(i - ROWS/2);
            x = round(j - COLS/2);
            x_cyl = s * atan(x/s); 
            y_cyl = (s * y) /sqrt(x^2 + s^2);
            
            % Move back to top left
            y_cyl = round(y_cyl + ROWS/2);
            x_cyl = round(x_cyl + COLS/2);
            
            % Assign Pixel values
            
            image_cylCoord(y_cyl, x_cyl, :) = image(i, j, :);
        end
    end

end