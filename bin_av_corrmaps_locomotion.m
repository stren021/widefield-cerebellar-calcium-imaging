
%averages, bins, and plots spatial maps for each locmotion epoch from the
%previous function

pixel_bin = -400:10:400;
pixel_bin = pixel_bin.*5;
pixel_bin = pixel_bin(1:(end-1));



for dom = 1
    
    within_all = COLLECT_MAPS_WITHIN;
    between_all = COLLECT_MAPS_BETWEEN;
    
    
    for wix = [1:4]
        
        within_curr = within_all(4:end, find(within_all(3,:) == wix));
        between_curr = between_all(4:end, find(between_all(3,:) == wix));
        
        within_av = nanmean(within_curr, 2);
        between_av = nanmean(between_curr,2);
        
        within_map = reshape(within_av, [801, 801]);
        between_map = reshape(between_av, [801, 801]);
        
            % Set bin (bin x bin square)
                bin = 10; % bin = 2 sets a 2x2 bin

                % Create empty matrix to store data in
                
                for pix = 1:2
                    switch pix
                        case 1
                            data = within_map;
                        case 2
                            data = between_map;
                    end
                
                       %binned_map = nan(size(data,1)/bin,size(data,2)/bin);
                       image_size = size(data,1);
                       bin_max = round(image_size/bin);


                        % will hold x position constant throughout all y pixels. acts as an
                        % anchor point for the bin
                        x_pos = 1;

                        for x = 1:bin_max
                            % Keeps x (row) constant for all y (column) pixels
                            x_pixel = x_pos;

                            % Set initial y (column) pixel
                            y_pos = 1;
                            for y = 1:bin_max
                                % Iteratively increases y pixel for bin based on bin size
                                y_pixel = y_pos;

                                % Averages pixels in bin (for 2 bin, will average a 2x2 square)
                                pixels_in_bin = nanmean(data(x_pixel:x_pixel+(bin-1),y_pixel:y_pixel+(bin-1)),'all');

                                % Assigns averaged value to corresponding binned pixel in new
                                % image
                                binned_image(x,y) = pixels_in_bin;

                                % Increases y pixel based on bin for next iteration
                                y_pos = y_pos+bin;
                            end

                            % Increases x pixel based on bin for next iteration
                            x_pos = x_pos+bin;
                        end

                        % Allocates new binned image to binned stack
                        binned_map = binned_image;
           
                        binned_map_x = nanmean(binned_map,1);
                        binned_map_y = nanmean(binned_map,2);
                        
                
                        if wix == 1

                            rest_map = binned_map;

                        else
                            diff_map = binned_map - rest_map;
         
                            imAlpha = ones(size(diff_map));
                            imAlpha(isnan(diff_map))=0;
                            figure(pix*10);subplot(2, 3, wix-1);imagesc(diff_map,'AlphaData',imAlpha);axis square;caxis([-0.15 0.15])
                            set(gca,'color',0.5*[1 1 1])
                            
                            

                            
                        end
                    
                end
                
                
                
                
                end 
                
                
    end
                
                
                
  