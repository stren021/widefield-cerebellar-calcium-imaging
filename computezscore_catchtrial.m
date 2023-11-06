%%%%computes z score maps aligned to auditory tone for the catch trial
dir_all = full_data_directory_here

file_name_out = [dir_all, 'zscores_beep_longwin'];

im_size = 1024;
mouse_list = dir(dir_all);
mouse_list = mouse_list(~ismember({mouse_list.name},{'.','..'}));

collect_peak = [];
collect_time = [];

fi = 1;
pos = 0;

bwin = 0:50:800;
tau_catch = 0:100:3000;

COLLECT_ALL_MAPS = [];




collect_onset = [];
collect_peak = [];

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

 compressed_name = ([dir_in, 'compressed_data.mat']);
 load (compressed_name);
% %%compressed data just has hU, hS, hV variables

events_name = ([dir_in, 'behavioral_events.mat']);
load (events_name);
%%this is the output mat from reaching_extract_triggers_v2

dend_som = ([dir_in, 'dend_som_IC_new.mat']);
load (dend_som);


 
% [imfileName,imfilePath] = uigetfile('*.tif','choose baseline image file'); % choose background image to plot ICs over
% im = imread([imfilePath,imfileName]);

%file_name_out = [dir_in, 'IC_dff_events_maps'];

cut_length = 10; %ten seconds removed from beginning of stack in my analysis
fs = 20; %20hz for imaging data
win = 2;
win_bin = -2000:50:2000;

total_length = (size(hV,1));
stack_num = size(beep_mat,2);

stack_length = total_length/stack_num;


event_mat = [];

IC_list = unique(IC_mask_mat(1,:));


for s_ix = 1:stack_num
    %this for loop just adjusts the timing of the events based on the amount cut off
    %per stack, and then makes an event_mat matrix with the timing indices
    %of the different events.
    
    event_mat_ix = [];
    bcurr = beep_mat(:,s_ix);
    
    rest = 1:400;
    
    if ~isnan(bcurr)
    
    bcurr = bcurr - cut_length;
    bcurr = round(bcurr.*fs);
    
    pmat = probmat(:,s_ix);
    p_omit = find(pmat == 2);
    pcurr = (bcurr(p_omit));
    
    event_mat_ix = zeros(stack_length,1);
    
    event_mat_ix(pcurr) = 1;
    event_mat_ix(rest) = -1;

    
    
    event_mat = [event_mat; event_mat_ix];
    end
 
 
end

    ctimes = find(event_mat == 1);
    rtimes = find(event_mat == -1);
    
    rtimes_2D = reshape(rtimes, [length(rest), length(rtimes)/length(rest)])';


    
    for dendsom = 1:2
        switch dendsom
            case 1
                IClist = dend;
            case 2
                IClist = som;
        end
        
        pos = 0;
        fi = fi + 1;
        
        DFF_MAP_CATCH = NaN(im_size^2,(length(tau_catch)-1));


        for IC_ix = 1:length(IClist)%(IC_list)1:length(som)%
            %pick your IC of interest here
            IC_curr = IClist(IC_ix);
            Reg_ID_2D = IC_mask_mat(4:end,IC_curr);

            Reg_id = reshape(Reg_ID_2D, [1024, 1024]);
            Reg_list = unique(Reg_ID_2D);
            Reg_list = Reg_list(2:end);

    
            DF_all = DF_ALL{IC_curr};
        
        
         

                     
         
         for dom_ix = 1:size(DF_all,1)
             
             dfcurr = DF_all(dom_ix,:);
             dom_ind = find(Reg_ID_2D == Reg_list(dom_ix));
             
             rest_mat = [];
             
             for r_ix = 1:size(rtimes_2D,1)
                 
                 rcurr = dfcurr(rtimes_2D(r_ix,:));
                 
                 rest_mat = [rest_mat; rcurr];
                 
             end
                 
             dfrest_av = mean(rest_mat,1);
             mean_rest = mean(dfrest_av);
             std_rest = std(dfrest_av);
             
             
             for event_type = 1
                 switch event_type
                     case 1
                         times = ctimes;
                         timeind = 1;
                         bwin = 0:50:3000;
                         taus = tau_catch;
                         dff_map = DFF_MAP_CATCH;
                 end
             



                        dfwin_collect = [];

                         for b_ix = 1:length(times);

                             bcurr = times(b_ix);

                             if event_type < 3
                                 beep_bin  = (bcurr):bcurr+(length(bwin)-1);
                             
                             else
                                 beep_bin = bcurr-8:bcurr+8;
                             end
                             
                             df_bin = dfcurr(beep_bin);

                             dfwin_collect = [dfwin_collect; df_bin];
                             
                             
                         end

                         dfwin_mean = mean(dfwin_collect);
                         z = (dfwin_mean-mean_rest)/std_rest;
                         
                         
                         for tau_ix = 1:(length(taus)-1)
                             start_bin = find(bwin == taus(tau_ix));
                             end_bin = find(bwin == taus(tau_ix + 1));
                             
                             zbin = z(start_bin:end_bin);
                             dfmean_bin = zbin(find(abs(zbin) == max(abs(zbin))));
                             
                             
                             
                             if abs(dfmean_bin) >= 2
                             dff_map(dom_ind,tau_ix) = dfmean_bin;
                             end
                             
                         end
                         
                 switch event_type
                     
                     case 1
                         DFF_MAP_CATCH = dff_map;

                 end
             end
         end
        end
       
        
        
    for pix = 1   
        switch pix
            case 1
                dff_map = DFF_MAP_CATCH;
                taus = tau_catch;
        end
        
        %fi = fi + 1;
        
        for tau_ix = 1:length(taus)-1
            
            mapcurr = dff_map(:,tau_ix);
            
            mapcurr_2D = reshape(mapcurr, [1024 1024]);
            
            COLLECT_ALL_MAPS = [COLLECT_ALL_MAPS [mouse_ix; day_ix; dendsom; pix; tau_ix; mapcurr]];
            
            imAlpha = ones(size(mapcurr_2D));
            imAlpha(isnan(mapcurr_2D)) = 0;
            
            pos = pos + 1;
            if pos > 16
                fi = fi + 1;
                pos = 1;
            end
            
            %figure(dendsom*10+pix);subplot(3,3,tau_ix);imshow(im, [1000 15000]);colormap gray;freezeColors;hold on;
            figure(fi);subplot(4,4,pos);imagesc(mapcurr_2D,'AlphaData',imAlpha);colormap parula;caxis([-3 3]);title(sprintf('t=%d',taus(tau_ix)))
            
        end
            
                
        
        
    end
    end
 
 end
end
end

save(file_name_out, 'COLLECT_ALL_MAPS', '-v7.3');