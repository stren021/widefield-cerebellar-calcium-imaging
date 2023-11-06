%%ICA only
%%takes a compressed_data.mat file and runs spatial ICA
%%file in structure should be mouseid/date/ e.g.
%%PcP2-041/08042022/compressed_data.mat




clear all

xDim = 1024; %image x-dimension (pixels)
yDim = 1024; %image y-dimension (pixels)

dir_all = 'F:\MSI\' %directory where all mouse data (compressed) is saved. so e.g. MSI\PcP2-041\etc
mouse_list = dir(dir_all);
mouse_list = mouse_list(~ismember({mouse_list.name},{'.','..'}));


for mouse_ix = 1 %mouse number in dir_all directory list (this is so you can batch if you need)
tic

mouse_dir = [mouse_list(mouse_ix).folder, '\' mouse_list(mouse_ix).name, '\']


day_list = dir(mouse_dir)
day_list = day_list(~ismember({day_list.name},{'.','..'}));

for day_ix = 1 %day number in mouse directory (this is so you can batch if you need)
   dir_in = [day_list(day_ix).folder, '\' day_list(day_ix).name, '\'] 



file_name_out = [dir_in 'ICA']

file_name = ([dir_in, 'compressed_data.mat']);

load (file_name)


    
    




B = jader_lsp([hU*hS]',60); %run spatial ICA using the jade algorithm 
%data must be in orientation where pixels are within a row and components
%are in each column (this is why you give it hU*hS'); second argument of
%the function is the number of ICs to look for; B is the unmixing matrix
%B = jader([hU*hS]',60);
S = hU*hS*B'; %multiply B into the original SVD data to get back your source signals (ICs)
for i = 1:size(S,2);
    sig = zscore(S(:,i));
    S_z(:,i) = sig;
end
ICs = reshape(S_z,[1024,1024,60]); %reshape your source signals into a 3D matrix of 256x256 images
% plot them using imagesc



%% Plot your ICs (note this is the raw IC z score map)
F = 1; %set figure name to 1
pos = 0; %set position in subplot to zero
for i = 1:size(ICs,3); %for each IC image
    pos = pos+1 %iterate position in variable pos by 1
    if pos > 16; %if position is > 16 set position to 1 and figure name to F+1;
        pos = 1;
        F = F+1;
        figure(F*1000);subplot(4,4,pos); imagesc(ICs(:,:,i)); pbaspect([1 1 1]); caxis([-3 3]); %plot ICs with -5 to 5 z-score on colorscale
    else
        figure(F*1000);subplot(4,4,pos); imagesc(ICs(:,:,i)); pbaspect([1 1 1]); caxis([-3 3]);
    end
end


save(file_name_out, 'B','S','ICs', '-v7.3');

end

end