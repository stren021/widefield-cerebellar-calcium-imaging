

%It takes an individual day's ICA run and thresholds all ICs
%It then plots a fig for each IC with the unthresholded and thresholded
%images
%Each figure number corresponds to IC number
%User goes through and identifies IC numbers that correspond to real ICs,
%excluding for vasculature etc.
%Put those IC numbers in a variable called IC_list_all, then sort them
%then re-run function with the file_name_out set to ICA_masks
%and ind_ix re-specified as IC_list_all


IC_mask_mat = [];
IC_mask_alldom = [];


dir_in = data_directory_here %directory where ICA.mat is stored

%%%%%%%%%%%%%%%%%%%%%%%%%%for all ICs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_name_out = [dir_in 'ICA_masks_all']


%%%%%%%%%for list of user identified ICs in variable: IC_list_all%%%%%%%%%%
%make sure you run command IC_list_all = sort(IC_list_all) prior ro running
%file_name_out = [dir_in 'ICA_masks']
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



file_name = ([dir_in, 'ICA']);

load (file_name)


pos = 0;

im_size = 1024;
mask=[1:im_size^2];
z_thr = 3.5;%3.5;
a_thr = 30;
mask=[1:im_size^2];
mapsig=zeros(size(mask));

IC_2D_all = zeros(im_size, im_size);

ic_flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%for all ICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ind_ix = 1:size(ICs,3)%
IC_ix = ind_ix;
%%%%%%%%%%%%%%%%%%%%%%%%for user identified ICs%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for ind_ix = 1:length(IC_list_all)%1:size(ICs,3)
%IC_ix = IC_list_all(ind_ix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig=S(:,IC_ix);
mapsig=sig;
map=reshape(mapsig,[im_size,im_size]);
map(abs(map)<z_thr)=0;

[Reg,FinSize,DomId] = ClustReg(abs(map),a_thr);

if ~isempty(DomId) %if DomId variable is populated
Reg_id=zeros(im_size,im_size); %create a blank image of zeros called Reg_id
%this part makes each domain in the image a different color
ic_flag = ic_flag + 1;
for ix=1:length(DomId) %for each entry in DomId indexed by ix
Reg0=Reg; % make Reg0 equal to Reg (from ClustReg)
id0=DomId(ix); %make id0 equal to the DomId value
Reg0(Reg0~=id0)=0; %make everything in Reg0 that's
%not equal to id0 equal to zero
Reg0(Reg0>0)=ix; %make anything greater than zero
%in Reg0 equal to ix (or equal to one value)
Reg_id=Reg_id+Reg0; %Sum the Reg0 image to Reg_id
%this makes each domain of an IC a different
%value/color




end

reg_reshape = reshape(Reg_id,[im_size^2],1);
reg_alldomain = zeros(size(reg_reshape));
ic_ind = find(reg_reshape > 0);
reg_alldomain(ic_ind) = ic_flag;

reg_alldomain_reshape = reshape(reg_alldomain, [im_size, im_size]);

IC_2D_all = IC_2D_all + reg_alldomain_reshape;

IC_mask_mat = [IC_mask_mat [IC_ix; length(DomId); mean(FinSize); reg_reshape]];
IC_mask_alldom = [IC_mask_alldom [IC_ix; reg_alldomain]];

end

figure(IC_ix*10);subplot(1,3,1);imagesc(Reg_id);pbaspect([1 1 1]); caxis([0 length(DomId)]);
figure(IC_ix*10);subplot(1,3,2);imagesc(ICs(:,:,IC_ix)); pbaspect([1 1 1]); caxis([-5 5]); %plot ICs with -5 to 5 z-score on colorscale
figure(IC_ix*10);subplot(1,3,3);imagesc(reg_alldomain_reshape);pbaspect([1 1 1]); caxis([0 ic_flag]);

end

figure(IC_ix*10 + 1); imagesc(IC_2D_all); pbaspect([1 1 1]);caxis([0 ic_flag]);

save(file_name_out, 'IC_mask_mat', 'IC_mask_alldom', 'IC_2D_all', '-v7.3');
