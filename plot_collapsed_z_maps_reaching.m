%%averages and plots the population z maps across all mice

tau_beep = 0:100:3000;
tau_drop = 0:100:800;
tau_touch = -400:100:400;

tau_mat = [tau_beep];

pos = 0;
fi = 20;


for dendsom = 1:2
    
    dendsom_ind = find(COLLECT_ALL_MAPS(3,:) == dendsom);
    
    dendsom_subset = COLLECT_ALL_MAPS(:,dendsom_ind);
    
    collect_av = NaN(length(tau_beep),1);
    collect_std = NaN(length(tau_beep),1);
    
        
        for tau_ix = 1:length(tau_beep)-1;
            
            tau_ind = find(dendsom_subset(5,:) == tau_ix);
            
            tau_subset = dendsom_subset(:,tau_ind);
            
            maps = tau_subset(6:end,:);
            
            map_av_2D = nanmean(maps,2);
            
            map_av_1D = nanmean(map_av_2D);
            map_std_1D = nanstd(map_av_2D);
            

            
            collect_av(tau_ix) = map_av_1D;
            collect_std(tau_ix) = map_std_1D;
            

        end
        
        figure(dendsom);hold on;shadedErrorBar(tau_beep, collect_av, collect_std, 'lineprops', 'y');
    end
