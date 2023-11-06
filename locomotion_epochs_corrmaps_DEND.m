%this progam creates domain centered maps of correlations with other
%domains either within a given IC or between a given IC. The window bin is
%in pixels, with 5 microns per pixel at the 1024x1024 resolution. A 401x401
%bin thus gives 2mm x 2mm

dir_all = full_data_directory_here
mouse_list = dir(dir_all);
mouse_list = mouse_list(~ismember({mouse_list.name},{'.','..'}));

file_name_out = [dir_all, 'corrmaps_locomotion_DEND_maps'];
COLLECT_AV_WITHIN = [];
COLLECT_AV_BETWEEN = [];
COLLECT_MAPS_WITHIN = [];
COLLECT_MAPS_BETWEEN = [];

stack_length_ca = 3400;
fs = 20;
winbin = 400; %window in pixels (100 = 500 microns) 
im_size = 1024;
stack_length = 13600;
stack_list = 4;
stack_id = 1:stack_list;
stack_id = stack_id.*stack_length_ca;
stack_id = [1 stack_id];

for mouse_ix = 1:length(mouse_list)
tic

mouse_dir = [mouse_list(mouse_ix).folder, '\' mouse_list(mouse_ix).name, '\']%'D:\MSI\PcP2-012\'


day_list = dir(mouse_dir)
day_list = day_list(~ismember({day_list.name},{'.','..'}));

for day_ix = 1:length(day_list) 
dir_in = [day_list(day_ix).folder, '\' day_list(day_ix).name, '\', 'stacks', '\', 'zoom1', '\']

file_name_segment = ([dir_in, 'locomotion_segmented_cat.mat']);

if exist(file_name_segment)
load(file_name_segment);

file_name = ([dir_in, 'ICA_masks.mat']);
load (file_name);

file_name_dff = ([dir_in, 'IC_dff.mat']);
load(file_name_dff);

file_name_dendsom = ([dir_in, 'dend_som_IC_new.mat']);
load(file_name_dendsom);



% COLLECT_MAPS_WITHIN = [];
% COLLECT_MAPS_BETWEEN = [];

win_walk = 3;
min_walk = 6;

locomat_all = zeros(stack_length,1);
locomat_rest = zeros(stack_length,1);
for stack_ix = 1:(length(stack_id)-1)
    
    walk_stack = start_stop_cat(find(start_stop_cat(:,1) > stack_id(stack_ix) & start_stop_cat(:,1) <= stack_id(stack_ix+1)),:);
    
    for walk_ix = 1:size(walk_stack,1);
    
        walk_onset_curr = walk_stack(walk_ix,1);
        walk_offset_curr = walk_stack(walk_ix,2);
        
        walk_bout = walk_onset_curr:walk_offset_curr;
        dur_walk_bout = length(walk_bout)/fs;
        
        locomat_rest(walk_bout) = 1;
        
        %%%%%%%takes only walk bouts > 5 sec%%%%%%%%%%%%
        if dur_walk_bout > min_walk
        
        %%%%%%%%%%%%%%%%%%pre walk window = 3 sec prior to walk onset%%%%%
        pre_walk = walk_onset_curr - win_walk*fs;
        
        %%%%%%%%%%%%%%%%%if 3 second window is longer than to the beginning
        %%%%%%%%%%%%%%%%%of the stack, it only takes to the begining of the
        %%%%%%%%%%%%%%%%%stack
        if pre_walk < stack_id(stack_ix)
            pre_walk = stack_id(stack_ix);
        end
        
        
        %%%make sure pre walk doesn't overlap with previous walk offset%%%
        %%%%%prev walk must be > 1 sec to avoid junk times%%%%%%%%%%%
        if walk_ix > 1
            
            walk_offset_prev = walk_stack(walk_ix-1,2);
            walk_onset_prev = walk_stack(walk_ix-1,1);
            walk_bout_prev = walk_onset_prev:walk_offset_prev;
            dur_walk_bout_prev = length(walk_bout_prev)/fs;
            
            if dur_walk_bout_prev > 1
                
            %%%walk offset plus 3 sec 'end walk' period%%%%%%%%%%%%%%%
            walk_offset_prev = walk_offset_prev + win_walk*fs;
            
                if pre_walk < walk_offset_prev
                    pre_walk = walk_offset_prev;
                end
                
            end
        end
        
        
        pre_walk_ind = pre_walk+1:walk_onset_curr;
    
        walk_bout = walk_onset_curr:walk_offset_curr;
        dur_walk_bout = length(walk_bout)/fs;
    
    
        %%%%%%%%%%%%defines initial waking as 3 sec window from walk onset
        if dur_walk_bout > win_walk
            initial_walk = walk_bout(1:win_walk*fs);
            
        %%%%%%%%%%calls continued walking everything beyond that 3 sec win    
            continued_walk = walk_bout(win_walk*fs+1:end);
            locomat_all(continued_walk) = 2;
            

        else
            initial_walk = walk_bout;
        end
        
        
        %%%%%%calls final walking the last 3 sec of walking%%%
        %%%%%note that this will overwrite any 'continued walking' window
        %%%%%overlap with final walking so that the continued walking is
        %%%%%independent of final walking
        
        if dur_walk_bout > win_walk*2
            final_walk = walk_bout((length(walk_bout)-((fs*win_walk)-1)):end);
            locomat_all(final_walk) = 3;
        end
        
        
        %%%%%%%%walk end is 3 seconds after the end of a walk bout%%%%%%%
        walk_end = walk_bout(end) + win_walk*fs;
        
        %%%%%%%%%%makes sure it doesn't overlap with the end of a stack%%
        if walk_end > stack_id(stack_ix + 1)
            walk_end = stack_id(stack_ix + 1);
        end
        
        %%%%%%makes sure walk end doesn't overlap with pre walk
        if walk_ix < size(walk_stack,1)
            walk_pre_next = walk_stack(walk_ix + 1,1) - win_walk*fs;
            if walk_end > walk_pre_next
                walk_end = walk_pre_next;
            end
        end
        
        
        post_walk_ind = (walk_bout(end)+1):walk_end;
        
        locomat_all(pre_walk_ind) = -1;
        locomat_rest(pre_walk_ind) = -1;
        locomat_all(initial_walk) = 1;
        locomat_all(post_walk_ind) = 4;
        locomat_rest(post_walk_ind) = 4;
    
    
    
    
    end

    end

end


pre_walktimes = find(locomat_all == -1);
initial_walktimes = find(locomat_all == 1);
continued_walktimes = find(locomat_all == 2);
final_walktimes = find(locomat_all == 3);
resttimes = find(locomat_rest == 0);
post_walktimes = find(locomat_all == 4);

length_tot = length(pre_walktimes) + length(initial_walktimes) + length(continued_walktimes) + length(final_walktimes) + length(post_walktimes) + length(resttimes)

fi = 1;
fi_b = 1;
pos = 0;
pos_b = 0

somind_between = 0;
dendind_between = 0;


IC_list = dend;
between_ind = somind_between;

map_all_between = [];



dff_collect = [];
index_collect = 0;
line_collect = [];
domain_centroid_index = [];

fi = fi*10;
fi_b = fi_b*10;

collect_x = [];
collect_y = [];

    for IC_ix = 1:length(IC_list); 

        IC_curr = DF_ALL{IC_list(IC_ix)};

        IC_map_2D = IC_mask_mat(4:end,IC_list(IC_ix));
        
        

        %centermat = zeros(im_size, im_size);
        
        for domain = 1:size(IC_curr,1)
            %this goes through each domain and finds the centroid for making
            %the correlaion maps
            
            domain_indices = find(IC_map_2D == domain);
            
            test = zeros(length(IC_map_2D),1);
            test(domain_indices) = domain;
            test_2D = reshape(test, [im_size im_size]);
            
            [x, y] = find(test_2D == domain);
            xdata=round(mean(y));
            ydata=round(abs(mean(x)));
            
            collect_x = [collect_x; xdata];
            collect_y = [collect_y; ydata];
            
            centermat = zeros(size(test_2D));
            
            centermat(ydata, xdata) = domain;
            
            centermat_1D = reshape(centermat, [1024^2,1]);
            domain_centroid = find(centermat_1D == domain);

            domain_centroid_index = [domain_centroid_index; domain_centroid];

        end

        %concatenate all DFF data to generate the large corellelogram
        dff_collect = [dff_collect; IC_curr];

        %make an index of the start rows of each new IC to keep track of what
        %domain belongs to what IC
        if IC_ix == 1
        index_collect = [index_collect; size(IC_curr,1)];
        line_collect = [line_collect; size(IC_curr,1)];
        else
        index_collect = [index_collect; (size(IC_curr,1)+index_collect(IC_ix))];
        line_collect = [line_collect; (size(IC_curr,1)+line_collect(IC_ix-1))];
        end
    end
    
    
    
    for wix = 1:6
        switch wix
            case 1
                segment = resttimes;
            case 2
                segment = pre_walktimes;
            case 3
                segment = initial_walktimes;
            case 4
                segment = continued_walktimes;
            case 5
                segment = final_walktimes;
            case 6
                segment = post_walktimes;
        end
        
        
        
        dff_select = dff_collect(:,segment);
    
        %make the correlation matrix for all domains across all ICs
        R_ALL = corrcoef(dff_select');

        %these are the structures that will hold each within/between IC correlation
        map_all_within = nan(((winbin)*2+1),((winbin)*2+1),size(R_ALL,1));

        map_all_between = nan(((winbin)*2+1),((winbin)*2+1),size(R_ALL,1));

        figure(mouse_ix*10+day_ix);subplot(3,2,wix);imagesc(R_ALL);axis square; pbaspect([1 1 1]); caxis([0 1])

            for l = 1:length(line_collect)-1

                l_ix = line_collect(l);
                figure(mouse_ix*10+day_ix);subplot(3,2,wix);hold on;line([l_ix, l_ix], [0, size(R_ALL,1)],'Color','k','LineWidth',0.75) 
                figure(mouse_ix*10+day_ix);subplot(3,2,wix);hold on;line([0, size(R_ALL,1)],[l_ix, l_ix],'Color','k','LineWidth',0.57) 

            end
            
            for IC_ix = 1:length(index_collect)-1

                %start and end indices for all the domains within an IC
                ind_start = index_collect(IC_ix)+1;
                ind_end = index_collect(IC_ix+1);
                ind_within = ind_start:ind_end;

                %start and end indices for all the domains between an IC
                ind_between = 1:max(index_collect);
                ind_between(ind_start:ind_end) = [];

                %select the parts of the correlelogram corresponding to all
                %within IC and between IC correlations for IC_ix. Columns are
                %always going to be for within the IC.
                R_within = R_ALL(ind_within,ind_within);
                R_between = R_ALL(ind_between,ind_within);
                

                collect_av_within = [];
                collect_av_between = [];

                for dom_ix = 1:size(R_within,2)

                    %makes blank 1D maps to put correlation values in
                    corrmap_within_2D = nan(im_size^2,1);
                    corrmap_between_2D = nan(im_size^2,1);

                    %finds the centroids for all domains within IC versus
                    %between IC for this particular domain
                    domain_centroid_within = domain_centroid_index(ind_within);
                    domain_centroid_between = domain_centroid_index(ind_between);

                    %finds corresponding correlation values
                    within_corrs = R_within(:,dom_ix);
                    between_corrs = R_between(:,dom_ix);
                    
                    collect_av_within = [collect_av_within; mean(within_corrs)];
                    collect_av_between = [collect_av_between; mean(between_corrs)];
                    

                    %puts the correlation values in the corresponding indices
                    %for the centroid of the domains
                    corrmap_within_2D(domain_centroid_within) = within_corrs;
                    corrmap_between_2D(domain_centroid_between) = between_corrs;

                    %centroid for the reference domain
                    dom_ix_centroid = domain_centroid_within(dom_ix);

                    %gets 2D coordinates for domain centroid
                    xy_centroid = zeros(im_size^2,1);
                    xy_centroid(dom_ix_centroid) = dom_ix;

                    xy_centroid = reshape(xy_centroid, [im_size im_size]);

                    [x, y] = find(xy_centroid == dom_ix);

                    xdata = round(mean(y));
                    ydata = round(abs(mean(x)));

                    %reshapes correlation maps into 2D structures
                    corrmap_within_2D = reshape(corrmap_within_2D, [im_size im_size]);
                    corrmap_between_2D = reshape(corrmap_between_2D, [im_size im_size]);


                    xmin = xdata-winbin;
                    xmax = xdata+winbin;
                    ymin = ydata-winbin;
                    ymax = ydata+winbin;

                    if xmin <= 0
                        xmin = 1;
                        xindmin = xdata-1;
                    else
                        xindmin = winbin;
                    end
                    if xmax > 1024
                        xmax = 1024;
                        xindmax = 1024-xdata;
                    else
                        xindmax = winbin;
                    end
                    if ymin <= 0
                        ymin = 1;
                        yindmin = ydata-1;
                    else
                        yindmin = winbin;
                    end
                    if ymax > 1024
                        ymax = 1024;
                        yindmax = 1024-ydata;
                    else
                        yindmax = winbin;
                    end

                    mapcurr_within = nan((winbin)*2+1,(winbin)*2+1);
                    mapcurr_between = nan((winbin)*2+1,(winbin)*2+1);

                    select_within = corrmap_within_2D(ymin:ymax,xmin:xmax);
                    select_between = corrmap_between_2D(ymin:ymax,xmin:xmax);

                    mapcurr_within((winbin+1)-yindmin:(winbin+1)+yindmax,(winbin+1)-xindmin:(winbin+1)+xindmax) = select_within;
                    mapcurr_between((winbin+1)-yindmin:(winbin+1)+yindmax,(winbin+1)-xindmin:(winbin+1)+xindmax) = select_between;

                    map_all_within(:,:,((ind_start-1)+dom_ix)) = mapcurr_within;
                    map_all_between(:,:,((ind_start-1)+dom_ix)) = mapcurr_between;

                end %dom_ix loop
                
            COLLECT_AV_WITHIN = [COLLECT_AV_WITHIN; [mouse_ix day_ix wix mean(collect_av_within)]];
            COLLECT_AV_BETWEEN = [COLLECT_AV_BETWEEN; [mouse_ix day_ix wix mean(collect_av_between)]];

            end %IC_ix loop


            mean_within = nanmean(map_all_within,3);


            mean_between = nanmean(map_all_between,3);
            
            
            mean_within_1D = reshape(mean_within, [((winbin)*2+1)^2,1]);
            mean_between_1D = reshape(mean_between, [((winbin)*2+1)^2,1]);
            
            COLLECT_MAPS_WITHIN = [COLLECT_MAPS_WITHIN [mouse_ix; day_ix; wix; mean_within_1D]];
            COLLECT_MAPS_BETWEEN = [COLLECT_MAPS_BETWEEN [mouse_ix; day_ix; wix; mean_between_1D]];
            

            
    end
end
end
end
    
save(file_name_out, 'COLLECT_MAPS_WITHIN', 'COLLECT_MAPS_BETWEEN');
