%% SVD_and_registration_test
%% this code is used to trouble shoot issues with lines showing up in our ICs 
% starts with a stack loading function and getting the indexes of blue and
% violet frames to de-interleave the stack

clear all



xDim = 1024; %image x-dimension (pixels)
yDim = 1024; %image y-dimension (pixels)
im_size = xDim*yDim;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%preprocessing params%%%%%%%%%%%%%%%%%%%%%%%%%%
stack_length = 13200;
split_length = stack_length\2;
cut_length = 10;
fs = 40;
fs_split = fs\2;
frame_cut = cut_length*fs_split;
corrected_length = split_length - frame_cut;
pcs2keep = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hData_cat= [];
mouse_dir = full_data_directory_here;
day_list = dir(mouse_dir)
day_list = day_list(~ismember({day_list.name},{'.','..'}));

for day_ix = 1:length(day_list)%2:length(day_list)
    
    day_dir = [day_list(day_ix).folder, '/' day_list(day_ix).name, '/']
    
    file_name_out = [day_dir 'compressed_data'];
    
    folder_dir=dir([day_dir, '0*']);


         hData_cat = zeros(im_size,(corrected_length*length(folder_dir)));

        for stack_ix = 1:length(folder_dir)
    
            tiff_dir = [folder_dir(stack_ix).folder, '/' folder_dir(stack_ix).name, '/']
            
            ind_start = ((stack_ix-1)*corrected_length) + 1;
            ind_end = ind_start + (corrected_length-1);

                %% Load in imaging stack and separate channels into blue and violet data matrices
                tic
                file_order=[];
                l=[];

                files_list=dir([tiff_dir, '*.tif']);



                for f=1:length(files_list)
                    fname=files_list(f).name;
                    l(f)=length(fname);
                end

                lmin=min(l);

                for f=1:length(files_list)
                    fname=files_list(f).name;

                    if length(fname)==lmin;
                        file_order(1)=f;
                    else
                        fseq=str2num(fname(end-8));
                        file_order(fseq+1)=f;
                    end
                end


    imname=fullfile(tiff_dir,files_list(1).name);
    im=tiffreadAltered(imname, [1], 'ReadUnknownTags',1);
    
    [yDim,xDim]=size(double(im.data))
    
    im_pixels=yDim*xDim
    
    data=nan(yDim,xDim,13200);
    
    fr=0;
  

                for f=file_order
                    fileName=fullfile(tiff_dir,files_list(f).name)

                    im_data=tiffreadAltered(fileName, [], 'ReadUnknownTags',1);

                    for cfr=1:length(im_data)
                        fr=fr+1;

                        data(:,:,fr)=double(im_data(cfr).data);
                    end
                end

                data(:,:,fr+1:end)=[];


            b_ind = [1:2:stack_length]; %get indexes of blue frames (for de-interleaving);
            v_ind = [2:2:stack_length]; %get indexes of violet frames (for de-interleaving);

            blue = data(:,:,b_ind); %get blue frame from stack variable using the b_ind variable of 
            % where the blue frames are put frames into 3D matrix blue where each page is an image


            blue = double(blue); % convert the 3D matrix of blue frames to a double
            blue = blue(:,:,frame_cut+1:split_length); %cut off the 1st 30sec of the file

            % for i = 1:length(v_ind); % for each violet frame 
            violet = data(:,:,v_ind); % get each violet image from the interleaved stack using the v_ind variable of 
            % %violet frame indices put these frames into a 3D matrix called violet where each page is an image
            % end
            % 
            violet = double(violet); % convert violet 3D matrix to a double
            violet = violet(:,:,frame_cut+1:split_length); % cut off the first 30 seconds of file


            bback = blue(:,:,round(split_length/2)); %get a blue image from the middle of the stack and call it bback
            fft_bback = fft2(bback); %perform a 2D fourier transform on the image and put into variable fft_bback

            vback = violet(:,:,round(split_length/2)); %get a violet image from the middle of the stack and call it vback
            fft_vback = fft2(vback); %perform a 2D fourier transform on the image and put into variable fft_vback

            clear data

            %%%%%%%%%%%%%%%%%%%%REGISTRATION AND FILTERING%%%%%%%%%%%%%%%%%%%%%%%%%%

            for t = 1:corrected_length; %for each frame of blue data (using variable t to index frames of blue stack)
            imb = blue(:,:,t); %get the frame and call it imb
            [BX(t,:),GregB] = dftregistration(fft_bback,fft2(imb),10); %register each frame by performing a fourier transform 
            %of the frame and  using the dftregistration function to register the image
            %to the previously calculated fourier transform fft_bback; last argument
            %in the function is the upsampling factor (this can be changed; 1 is no
            %upsampling)
            imb = abs(ifft2(GregB)); % take the absolute value of the inverse fourier transform
            %of the registered blue image
            %imb = imgaussfilt(imb,1,'FilterDomain','spatial'); % apply gaussian filter to the image
            blue_align(:,:,t) = imb; % put newly aligned, filtered images into to a new 3D matrix called blue_align
            end

            for t = 1:corrected_length; %for each frame of violet data (using variable t to index frames of violet stack)
            imv = violet(:,:,t); %get the frame and call it imv
            [VX(t,:),GregV] = dftregistration(fft_vback,fft2(imv),10); % register each frame by performing a fourier transform
            %of the frame and using the dftregistration function to register the image
            %to the previously calculated fourier transform fft_vback; last argument is
            %upsample factor as before with blue data
            imv = abs(ifft2(GregV)); % take the absolute value of the inverse fourier transform
            %of the registered violet image
            %imv = imgaussfilt(imv,1,'FilterDomain','spatial'); % apply gaussian filter to image
            violet_align(:,:,t) = imv; %put newly aligned,filtered images into new 3D matrix called violet_align
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%ALIGNMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            blue_2D_align = reshape(blue_align,[im_size,corrected_length]); %reshape aligned blue image stack to 2-dimensions
            %now each column of the image is a reshaped 1-D slice of the stack and time
            %is running in the columns
            violet_2D_align = reshape(violet_align,[im_size,corrected_length]);%reshape aligned violet image stack to 2-dimensions
            %now each column of the image is a reshaped 1-D slice of the stack and time
            %is running in the columns

            %%%%%%%%%%%%%%%%%%%%%%%%SVD ON BLUE AND VIOLET%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  [bU,bS,bV] = compute_svd(blue_2D_align,'randomized',200); %compute the singular values of the blue and violet data
            % %U is space, S is a diagonal matrix of singlular values, and V is time,
            % %randomized is the method of decomposition and 200 is the number of
            % %components kept by the SVD (rank of the data);
            %  [vU,vS,vV] = compute_svd(violet_2D_align,'randomized',200); 
            %  save('Cmprs_stk003_splitalign10_filt1spatial_030421','vU','vS','vV','bU','bS','bV','GregB','GregV','BX','VX');
            % % save the SVD of the two image stacks
            %% hemo correction
            
            
             bData= blue_2D_align;%bU*bS*bV'; %reconstruct the 2D blue and violet data by multiplying out the SVD
             clear blue_2D_align;
             clear blue_align
             clear blue

             vData= violet_2D_align;%vU*vS*vV';
             clear violet_2D_align;
             clear violet_align
             clear violet


             bData2 = detrend(bData'); 
             vData2 = detrend(vData');
             bData2 = bData2'; 
             vData2 =vData2';
            % 
             avb = mean(bData,2); % get the average blue and violet images over time (this should come out as a single image of 65536 pixels)
             avv = mean(vData,2);
            % 
             hData = bData2./avb-vData2./avv; %hemo correct the data by subtracting the detrended vData 
            % %divided by the average vData over time from the detrended bData divided by
            % %the average bData over time and put it into the variable hData
            
            hData_cat(1:im_size, ind_start:ind_end) = hData;
           
            clear hData
            clear bData2
            clear vData2
            clear vData
            clear bData
            clear data

        end
        
            [hU,hS,hV] = compute_svd(hData_cat,'randomized',pcs2keep); % SVD the hemo corrected data for easy loading later

            save(file_name_out, 'hU','hS','hV', '-v7.3');


            clear hU
            clear hS
            clear hV
    end

