%this code takes the compressed data and IC masks and generates the time
%series df/f for each IC

dir_in = put_data_location_here;


file_name = ([dir_in, 'ICA_masks.mat']);

 load (file_name);

file_name_compressed = ([dir_in, 'compressed_data.mat']);


 load(file_name_compressed);

 
file_name_out = [dir_in, 'IC_dff'];

 hData = hU*hS*hV';


IC_list = unique(IC_mask_mat(1,:));

DF_ALL = cell(1,length(IC_list));

for c_ix = 1:length(IC_list);
    DF_ALL(c_ix) = {[]};
end

for ix = 1:length(IC_list)

IC_ix = IC_list(ix);

Reg_ID_2D = IC_mask_mat(4:end,ix);


Reg_id = reshape(Reg_ID_2D, [1024, 1024]);
Reg_list = unique(Reg_ID_2D);
Reg_list = Reg_list(2:end);





DF_all = zeros(length(Reg_list),size(hV,1));

for dom_ix = 1:length(Reg_list)
    
    dom_ind = find(Reg_ID_2D == Reg_list(dom_ix));
    
    df_ix = hData(dom_ind,:);
    
    df_mean = mean(df_ix);
    
    DF_all(dom_ix,:) = df_mean;
    
end

timepoints = 1:size(DF_all,2);
timepoints = timepoints./20;

%DF_all = DF_all';

timepoints = 1:size(DF_all,2);
timepoints = timepoints./20;

for p = 1:size(DF_all,1);
dfcurr = DF_all(p,:);
dfcurr = dfcurr+((p-1)*0.2);
figure(IC_ix);hold on;plot(timepoints,dfcurr);xlim([0 170]);
end

DF_ALL(ix) = {[DF_all]};

figure(IC_ix*100); imagesc(Reg_id);pbaspect([1 1 1]); caxis([0 length(Reg_list)])
%figure(IC_ix*100+1);hold on; subplot(1,2,2);imagesc(ICs(:,:,IC_ix)); pbaspect([1 1 1]); caxis([-3 3]); %plot 
end

save(file_name_out, 'DF_ALL', '-v7.3');
