
%%%%collects average df/f for each locomotion segment%%%%%%%%%%

collect_df_som = [];
collect_df_dend = [];

dir_all = full_data_directory_here;
mouse_list = dir(dir_all);
mouse_list = mouse_list(~ismember({mouse_list.name},{'.','..'}));
%directory with all of the data


stack_length_ca = 3400;
fs = 20;
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




win_walk = 3;

%rest = 0
%pre walk = -1
%walk = 1
%continued walk = 2

locomat_all = zeros(stack_length,1);


for stack_ix = 1:(length(stack_id)-1)
    
    walk_stack = start_stop_cat(find(start_stop_cat(:,1) > stack_id(stack_ix) & start_stop_cat(:,1) <= stack_id(stack_ix+1)),:);
    
    for walk_ix = 1:size(walk_stack,1);
    
        walk_onset_curr = walk_stack(walk_ix,1);
        walk_offset_curr = walk_stack(walk_ix,2);
        
        walk_bout = walk_onset_curr:walk_offset_curr;
        dur_walk_bout = length(walk_bout)/fs;
        
        if dur_walk_bout > 3
        
            pre_walk = walk_onset_curr - win_walk*fs;

            if pre_walk < stack_id(stack_ix)
                pre_walk = stack_id(stack_ix);
            end

            if walk_ix > 1

                walk_offset_prev = walk_stack(walk_ix-1,2);

                if pre_walk < walk_offset_prev
                    pre_walk = walk_offset_prev;
                end
            end

            pre_walk_ind = pre_walk:walk_onset_curr;

            walk_bout = walk_onset_curr:walk_offset_curr;
            dur_walk_bout = length(walk_bout)/fs;



            if dur_walk_bout > win_walk
                initial_walk = walk_bout(1:win_walk*fs);
                continued_walk = walk_bout(win_walk*fs+1:end);

                locomat_all(continued_walk) = 2;

            else
                initial_walk = walk_bout;
            end
            
        if dur_walk_bout > win_walk*2
            final_walk = walk_bout((length(walk_bout)-((fs*win_walk)-1)):end);
            locomat_all(final_walk) = 3;
        end
            walk_end = walk_bout(end) + win_walk*fs;

            if walk_end > stack_id(stack_ix + 1)
                walk_end = stack_id(stack_ix + 1);
            end

            if walk_ix < size(walk_stack,1)
                walk_pre_next = walk_stack(walk_ix + 1,1) + win_walk*fs;
                if walk_end > walk_pre_next
                    walk_end = walk_pre_next
                end
            end


            post_walk_ind = walk_bout(end):walk_end;

            locomat_all(pre_walk_ind) = -1;
            locomat_all(initial_walk) = 1;
            locomat_all(post_walk_ind) = 4;
            
        end
    
    
    
    
    end



end

pre_walktimes = find(locomat_all == -1);
initial_walktimes = find(locomat_all == 1);
continued_walktimes = find(locomat_all == 2);
final_walktimes = find(locomat_all == 3);
resttimes = find(locomat_all == 0);
post_walktimes = find(locomat_all == 4);

fi = 1;
fi_b = 1;
pos = 0;
pos_b = 0

somind_between = 0;
dendind_between = 0;

for pix = 1:2
    
    switch pix
        case 1
            IC_list = som;
            collect_df = collect_df_som;
        case 2
            IC_list = dend;
            collect_df = collect_df_dend;
    end
    



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
        Reg_list = unique(IC_map_2D);
        Reg_list = Reg_list(2:end);
        
    
    
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
        

        dff_select = IC_curr(:,segment);
        df_mn = mean(dff_select(:,:));
        df_mn = mean(df_mn);
        
        collect_df = [collect_df; [wix df_mn]];
    

            
        
    
    end
     
    end
    

    
         switch pix
            case 1
                collect_df_som = collect_df;

            case 2
                collect_df_dend = collect_df;
         end
end
end
end
end %end PIX loop