%sliding window correlation analysis for reach vs no reach in standard
%paradigm


dir_all = full_data_directory_here

file_name_out = [dir_all, 'corrmovav_reachnoreach_SOM_2sec'];

im_size = 1024;
mouse_list = dir(dir_all);
mouse_list = mouse_list(~ismember({mouse_list.name},{'.','..'}));


cut_length = 10; %ten seconds removed from beginning of stack in my analysis
fs = 20;
pre_beep = 0.300;
post_beep = 3.7;

pre_touch = 3;
post_touch = 1;

pre_touch_win = pre_touch*fs;
post_touch_win = post_touch*fs;


pre_beep_win = pre_beep*fs;
post_beep_win = post_beep*fs;

reach_taus = -0.3:0.05:3.7;
touch_taus = -3:0.05:1;

slide_taus = -0.3:0.3:3.7;
plot_taus = -0.3:0.05:3.45;
bin_size = 6;



collect_av_within_reach = [];
collect_av_between_reach = [];
collect_av_within_noreach = [];
collect_av_between_noreach = [];

collect_z_within_reach = [];
collect_z_between_reach = [];
collect_z_within_noreach = [];
collect_z_between_noreach = [];

COLLECT_ALL_WITHIN_REACH = []
COLLECT_ALL_BETWEEN_REACH = [];
COLLECT_ALL_WITHIN_NOREACH = [];
COLLECT_ALL_BETWEEN_NOREACH = [];



for mouse_ix = 1:length(mouse_list)
tic

mouse_dir = [mouse_list(mouse_ix).folder, '\' mouse_list(mouse_ix).name, '\']


day_list = dir(mouse_dir)
day_list = day_list(~ismember({day_list.name},{'.','..'}));

for day_ix = 1:length(day_list) 
    
dir_in = [day_list(day_ix).folder, '\' day_list(day_ix).name, '\']


 dff_name = ([dir_in, 'IC_dff.mat']);
 

 
 if exist(dff_name) 
 load (dff_name);
 
mask_name = ([dir_in, 'ICA_masks.mat']);
load (mask_name);
%%ICA_masks.mat is a matrix with each column as an ICA mask, where each
%%domain has an integer index. The first 3 rows are the IC #, the # of
%%domains, and the average size of each domain. I remove those first three
%%rows below (line 32) so you can modify to have it work with a different
%%matrix format

collect_within_reach = [];
collect_between_reach = [];
collect_within_noreach = [];
collect_between_noreach = [];
% %%compressed data just has hU, hS, hV variables

events_name = ([dir_in, 'behavioral_events.mat']);
load (events_name);
%%this is the output mat from reaching_extract_triggers_v2

dend_som = ([dir_in, 'dend_som_IC_new.mat']);
load (dend_som);

 compressed_name = ([dir_in, 'compressed_data.mat']);
 load (compressed_name);
% %%compressed data just has hU, hS, hV variables

total_length = (size(hV,1));
stack_num = size(beep_mat,2);

stack_length = total_length/stack_num;

beep_window_reach = [];
beep_window_noreach = [];
rest_window_all = [];

touch_window = [];
latency2touch = [];


collect_av_within_rest = [];
collect_av_between_rest = [];

for s_ix = 1:stack_num
    %this for loop just adjusts the timing of the events based on the amount cut off
    %per stack, and then makes an event_mat matrix with the timing indices
    %of the different events.
    

    
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
    rest_window_all = [rest_window_all rest_window];
    
    for beep_ix = 1:length(bcurr)
        
        if beep_ix < length(bcurr)
        
            t_ind = find(tcurr > dcurr(beep_ix) & tcurr < bcurr(beep_ix + 1));
            
            touch_ind = tcurr(t_ind);
        else
            
            t_ind = find(tcurr > dcurr(beep_ix) & tcurr < stack_length);
            touch_ind = tcurr(t_ind);
            
        end
            
            if ~isempty(touch_ind)
                
                reach_win = (dcurr(beep_ix):touch_ind(1))';
                
                latency2touch = [latency2touch; length(reach_win)/fs];
                
                if (length(reach_win)/fs) <= 2
                    
                    b_window = bcurr(beep_ix)-pre_beep_win:bcurr(beep_ix)+post_beep_win;
                    stack_correct = (s_ix-1)*stack_length;
                    b_window = b_window + stack_correct;
                    
                    beep_window_reach = [beep_window_reach; b_window];
                    
                end
            else
                
                b_window = bcurr(beep_ix)-pre_beep_win:bcurr(beep_ix)+post_beep_win;
                stack_correct = (s_ix-1)*stack_length;
                b_window = b_window + stack_correct;
                beep_window_noreach = [beep_window_noreach; b_window];
                
            end
            
    end
end
            
            
            IC_list = som;
            dff_collect = [];
            line_collect = [];
            index_collect = 0;
            
            %%%%%%determine indices for IC start end for
            %%%%%%within/between
            
           for IC_ix = 1:length(IC_list);
                
                IC = IC_list(IC_ix);
                IC_map_2D = IC_mask_mat(4:end,IC);

                Reg_id = reshape(IC_map_2D, [1024, 1024]);
                Reg_list = unique(IC_map_2D);
                Reg_list = Reg_list(2:end);
                
                IC_curr = DF_ALL{IC};
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
            
           %%%%%rest correlation calculations%%%%%
                dff_select_rest = dff_collect(:,rest_window_all);
                R_ALL_REST = corrcoef(dff_select_rest');
                collect_rest_within = [];
                collect_rest_between = [];
                
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
                        R_within = R_ALL_REST(ind_within,ind_within);
                        R_between = R_ALL_REST(ind_between,ind_within);
                        
                        for dom_ix = 1:size(R_within,2)
                            
                            %finds corresponding correlation values
                            within_corrs = R_within(:,dom_ix);
                            between_corrs = R_between(:,dom_ix);
                        
                            collect_rest_within = [collect_rest_within; mean(within_corrs)];
                            collect_rest_between = [collect_rest_between; mean(between_corrs)];
                        end
                        
                    end
                     
                    collect_av_within_rest = [collect_av_within_rest; (collect_rest_within)];
                    collect_av_between_rest = [collect_av_between_rest; (collect_rest_between)];
           
                %%%%%%end rest calculations%%%%%%
                for pix = 1:2
                    switch pix
                        case 1
                            time_mat = beep_window_reach;
                            taus = reach_taus;
                            collect_all_within = collect_within_reach;
                            collect_all_between = collect_between_reach;
                            
                        case 2
                            time_mat = beep_window_noreach;
                            taus = reach_taus;
                            collect_all_within = collect_within_noreach;
                            collect_all_between = collect_between_noreach;

                 end
                    
                     collect_within = [];
                     collect_between = [];
                     
                     pos = 0;
                     

                     
                     for tau_ix = 1:(size(time_mat,2)-(bin_size-1))
                        
                        time_ix = time_mat(:,(tau_ix:(tau_ix+(bin_size-1))))';
                        time_ix_1D = reshape(time_ix, [1,(size(time_ix,1)*(size(time_ix,2)))]);
                        
                        dff_select = dff_collect(:,time_ix_1D);
                        R_ALL = corrcoef(dff_select');
                        
%                         pos = pos + 1;
%                     figure(pix);subplot(4,4,pos);imagesc(R_ALL);axis square; pbaspect([1 1 1]); caxis([0 1])
% 
%                     for l = 1:length(line_collect)-1
% 
%                         l_ix = line_collect(l);
%                         figure(pix);subplot(4,4,pos);hold on;line([l_ix, l_ix], [0, size(R_ALL,1)],'Color','k','LineWidth',1.5) 
%                         figure(pix);subplot(4,4,pos);hold on;line([0, size(R_ALL,1)],[l_ix, l_ix],'Color','k','LineWidth',1.5) 
% 
%                     end  
                    
                    
                        collect_av_within = [];
                        collect_av_between = [];
                    
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
                        
                        for dom_ix = 1:size(R_within,2)
                            
                            %finds corresponding correlation values
                            within_corrs = R_within(:,dom_ix);
                            between_corrs = R_between(:,dom_ix);
                        
                            collect_av_within = [collect_av_within; mean(within_corrs)];
                            collect_av_between = [collect_av_between; mean(between_corrs)];
                        end
                        
                     end
                    
                    collect_within = [collect_within (collect_av_within)];
                    collect_between = [collect_between (collect_av_between)];
                    
                    end
                    
                    collect_all_within = [collect_all_within; collect_within];
                    collect_all_between = [collect_all_between; collect_between];
                    
                switch pix
                    
                    case 1
                        collect_within_reach = collect_all_within;
                        collect_between_reach = collect_all_between;
                    case 2
                        collect_within_noreach = collect_all_within;
                        collect_between_noreach = collect_all_between;
                        
                end
                        
                end
                
                
 mean_rest_within = mean(collect_av_within_rest);
 std_rest_within = std(collect_av_within_rest);
 
 mean_rest_between = mean(collect_av_between_rest);
 std_rest_between = std(collect_av_between_rest);
 
 within_av_reach = mean(collect_within_reach);
 between_av_reach = mean(collect_between_reach);
 
 within_av_noreach = mean(collect_within_noreach);
 between_av_noreach = mean(collect_between_noreach);
 
 collect_av_within_reach = [collect_av_within_reach [mouse_ix; day_ix; within_av_reach']];
 collect_av_between_reach = [collect_av_between_reach [mouse_ix; day_ix; between_av_reach']];
 
 if ~isempty(beep_window_noreach);
 collect_av_within_noreach = [collect_av_within_noreach [mouse_ix; day_ix; within_av_noreach']];
 collect_av_between_noreach = [collect_av_between_noreach [mouse_ix; day_ix; between_av_noreach']];
 end
 
 within_z_reach = (collect_within_reach./collect_av_within_rest);
 between_z_reach = (collect_between_reach./collect_av_between_rest);
 within_z_reach_mean = mean(collect_within_reach./collect_av_within_rest);
 between_z_reach_mean = mean(collect_between_reach./collect_av_between_rest);
 
 if ~isempty(beep_window_noreach);
 within_z_noreach = (collect_within_noreach./collect_av_within_rest);
 between_z_noreach = (collect_between_noreach./collect_av_between_rest);
 within_z_noreach_mean = mean(collect_within_noreach./collect_av_within_rest);
 between_z_noreach_mean = mean(collect_between_noreach./collect_av_between_rest);
 end
 

collect_z_within_reach = [collect_z_within_reach [mouse_ix; day_ix; within_z_reach_mean']];
collect_z_between_reach = [collect_z_between_reach [mouse_ix; day_ix; between_z_reach_mean']];

if ~isempty(beep_window_noreach);
collect_z_within_noreach = [collect_z_within_noreach [mouse_ix; day_ix; within_z_noreach_mean']];
collect_z_between_noreach = [collect_z_between_noreach [mouse_ix; day_ix; between_z_noreach_mean']];
end


 COLLECT_ALL_WITHIN_REACH = [COLLECT_ALL_WITHIN_REACH; within_z_reach];
 COLLECT_ALL_BETWEEN_REACH = [COLLECT_ALL_BETWEEN_REACH; between_z_reach];
 
 if ~isempty(beep_window_noreach);
 COLLECT_ALL_WITHIN_NOREACH = [COLLECT_ALL_WITHIN_NOREACH; within_z_noreach];
 COLLECT_ALL_BETWEEN_NOREACH = [COLLECT_ALL_BETWEEN_NOREACH; between_z_noreach];
 end
                
 end
 

            
            
end
end

save(file_name_out, 'COLLECT_ALL_WITHIN_REACH', 'COLLECT_ALL_BETWEEN_REACH', 'collect_z_within_reach',...
    'collect_z_between_reach', 'collect_av_within_reach','collect_av_between_reach',...
    'COLLECT_ALL_WITHIN_NOREACH', 'COLLECT_ALL_BETWEEN_NOREACH', 'collect_z_within_noreach',...
    'collect_z_between_noreach', 'collect_av_within_noreach','collect_av_between_noreach');
            