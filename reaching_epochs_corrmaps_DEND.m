%%%%segments data and creates correlation maps similar to locomotion
%%%%segement analysis


dir_all = full_data_directory_here

im_size = 1024;
mouse_list = dir(dir_all);
mouse_list = mouse_list(~ismember({mouse_list.name},{'.','..'}));

file_name_out = [dir_all, 'corrmaps_reaching_DEND_maps'];

collect_peak = [];
collect_time = [];

fi = 1;
pos = 0;

bwin = 0:50:800;

collect_onset = [];
collect_peak = [];

latency2touch = [];
winbin = 400;
COLLECT_AV_WITHIN = [];
COLLECT_AV_BETWEEN = [];
COLLECT_MAPS_WITHIN = [];
COLLECT_MAPS_BETWEEN = [];


for mouse_ix = 1:length(mouse_list)
tic

mouse_dir = [mouse_list(mouse_ix).folder, '\' mouse_list(mouse_ix).name, '\']'


day_list = dir(mouse_dir)
day_list = day_list(~ismember({day_list.name},{'.','..'}));

for day_ix = 1:length(day_list) 
    
dir_in = [day_list(day_ix).folder, '\' day_list(day_ix).name, '\']



 dff_name = ([dir_in, 'IC_dff.mat']);
 
 if exist(dff_name)
     
 load (dff_name);


events_name = ([dir_in, 'behavioral_events.mat']);
load (events_name);
%%this is the output mat from reaching_extract_triggers_v2


mask_name = ([dir_in, 'ICA_masks.mat']);
load (mask_name);
%%ICA_masks.mat is a matrix with each column as an ICA mask, where each
%%domain has an integer index. The first 3 rows are the IC #, the # of
%%domains, and the average size of each domain. I remove those first three
%%rows below (line 32) so you can modify to have it work with a different
%%matrix format

 compressed_name = ([dir_in, 'compressed_data.mat']);
 load (compressed_name);
% %%compressed data just has hU, hS, hV variables



dend_som = ([dir_in, 'dend_som_IC_new.mat']);
load (dend_som);


%file_name_out = [dir_in, 'IC_dff_events_maps'];

cut_length = 10; %ten seconds removed from beginning of stack in my analysis
fs = 20; %20hz for imaging data
win = 2;
win_bin = -2000:50:2000;
pre_win = 5;
pre_win_bin = pre_win*fs;
post_win = 2;
post_win_bin = post_win*fs;

post_reward_dur = 5;
post_reward_bin = post_reward_dur*fs;


total_length = (size(hV,1));
stack_num = size(beep_mat,2);

stack_length = total_length/stack_num;


event_mat = [];

IC_list = unique(IC_mask_mat(1,:));

pre_reach = [];
pre_noreach = [];
cue_reach = [];
cue_noreach = [];
drop_reach = [];
drop_noreach = [];
reward = [];
post_reward = [];
rest = [];


noreach = 0;


for s_ix = 1:stack_num
    %this for loop just adjusts the timing of the events based on the amount cut off
    %per stack, and then makes an event_mat matrix with the timing indices
    %of the different events.
    
    event_mat_ix = [];
    
    
    bcurr = beep_mat(:,s_ix);
    bcurr = bcurr - cut_length;
    bcurr = round(bcurr.*fs);
    
    dcurr = drop_mat(:,s_ix);
    dcurr = dcurr - cut_length;
    dcurr = round(dcurr.*fs);
    
    tcurr = touch_mat(:,s_ix);
    tcurr = tcurr(find(tcurr > 0));
    tcurr = tcurr - cut_length;
    tcurr = round(tcurr.*fs);
        
    rest_window = 1:400;
    stack_correct = (s_ix-1)*stack_length;
    rest_window = rest_window + stack_correct;
    rest = [rest; rest_window];
    
     for beep_ix = 1:length(bcurr)
        
        beep_win = (bcurr(beep_ix):dcurr(beep_ix))';
        stack_correct = (s_ix-1)*stack_length;
        beep_win = beep_win + stack_correct;
        pre_win = ((bcurr(beep_ix)-pre_win_bin):bcurr(beep_ix))';
        pre_win = pre_win + stack_correct;
        
        if beep_ix < length(bcurr)
        
            t_ind = find(tcurr > dcurr(beep_ix) & tcurr < bcurr(beep_ix + 1));
            
            touch_ind = tcurr(t_ind);
        else
            
            t_ind = find(tcurr > dcurr(beep_ix) & tcurr < stack_length);
            touch_ind = tcurr(t_ind);
            
        end
            
            if ~isempty(touch_ind)
                
                reach_win = (dcurr(beep_ix):touch_ind(1))';
                reach_win = reach_win + stack_correct;
                
                latency2touch = [latency2touch; length(reach_win)/fs];
                
                if (length(reach_win)/fs) <= 2
                

                
                
                pre_reach = [pre_reach; pre_win];
                cue_reach = [cue_reach; beep_win];
                drop_reach = [drop_reach; reach_win];

                
                end
                
            else
                noreach = noreach + 1;
                
                pre_noreach = [pre_noreach; pre_win];
                cue_noreach = [cue_noreach; beep_win];

                noreach_win = (dcurr(beep_ix):dcurr(beep_ix)+post_win_bin)';
                noreach_win = noreach_win + stack_correct;
                drop_noreach = [drop_noreach; noreach_win];
                
            end
            
    end
    
end


if noreach > 0
    
    event_mat = [1 2 3 4 6 7 8];
else 
    event_mat = [1 2 3 4];
    
end
                

somind_between = 0;
dendind_between = 0;


        IC_list = dend;
        between_ind = somind_between;

        
        pos = 0
        
        map_all_between = [];
        dff_collect = [];
        index_collect = 0;
        line_collect = [];
        domain_centroid_index = [];
        


        collect_x = [];
        collect_y = [];
        
        
            for IC_ix = 1:length(IC_list); 

            IC = IC_list(IC_ix);
            IC_map_2D = IC_mask_mat(4:end,IC);

            Reg_id = reshape(IC_map_2D, [1024, 1024]);
            Reg_list = unique(IC_map_2D);
            Reg_list = Reg_list(2:end);
        
            IC_curr = DF_ALL{IC};
            

        
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
    
             for wix = event_mat;
                 switch wix
                     case 1
                         segment = rest;
                     case 2
                         segment = pre_reach;
                     case 3
                         segment = cue_reach;
                     case 4
                         segment = drop_reach;
                     case 6
                         segment = pre_noreach;
                     case 7
                         segment = cue_noreach;
                     case 8
                         segment = drop_noreach;
                 end
                 

        

        dff_select = dff_collect(:,segment);
        
                %make the correlation matrix for all domains across all ICs
        R_ALL = corrcoef(dff_select');

        %these are the structures that will hold each within/between IC correlation
        map_all_within = nan(((winbin)*2+1),((winbin)*2+1),size(R_ALL,1));

        map_all_between = nan(((winbin)*2+1),((winbin)*2+1),size(R_ALL,1));

        figure(mouse_ix*10+day_ix);subplot(4,2,wix);imagesc(R_ALL);axis square; pbaspect([1 1 1]); caxis([0 1])

            for l = 1:length(line_collect)-1

                l_ix = line_collect(l);
                figure(mouse_ix*10+day_ix);subplot(4,2,wix);hold on;line([l_ix, l_ix], [0, size(R_ALL,1)],'Color','k','LineWidth',0.75) 
                figure(mouse_ix*10+day_ix);subplot(4,2,wix);hold on;line([0, size(R_ALL,1)],[l_ix, l_ix],'Color','k','LineWidth',0.75) 

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

                    %rotates clockwise 90 degrees (something during reshaping
                    %process makes the ending 2D image rotated 90 degrees
                    %counter clockwise, I don't know why
                    %corrmap_within_2D = rot90(corrmap_within_2D,3);
                    %corrmap_between_2D = rot90(corrmap_between_2D,3);


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
end %day_ix loop
 %mouse_ix loop

save(file_name_out, 'COLLECT_MAPS_WITHIN', 'COLLECT_MAPS_BETWEEN', 'COLLECT_AV_WITHIN', 'COLLECT_AV_BETWEEN');