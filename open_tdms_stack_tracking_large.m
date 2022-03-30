% function [stack]=open_tdms_stack(save_name);

path='E:\current_data\M1-M2 project\20201022_m1234_m2mycotest\shape bf\m2_5\';
save_name='all_frames.tif';

files=dir([path,'*.tdms']);

% suppose 's' is the struct array. 'DOB' is the field that contains date and time.
T = struct2table(files); % convert the struct array to a table
sortedT = sortrows(T, 'date'); % sort the table by 'DOB'
files = table2struct(sortedT); % change it back to struct array if necessary
%%
Nframes=length(files);
Nfr2save=Nframes-1;
initial=Nframes-Nfr2save;

I1=get_image_from_tdms([path,files(initial).name]);
I1=flipud(rot90(I1));
I1=uint16(I1);

%%
imwrite(I1,[path,save_name]);
%
Nrows=size(I1,1);
Ncolumns=size(I1,2);
stack=zeros(Nrows,Ncolumns,Nfr2save);
stack(:,:,1)=I1;

%%
    ii=2;
for k=initial+1:Nframes

    I=get_image_from_tdms([path,files(ii).name]);
    I=flipud(rot90(I));
%     I=imcrop(I,r);
    I=uint16(I);
    stack(:,:,ii)=I;
    imwrite(I,[path,save_name],'WriteMode','append');
    disp(['done ',num2str(k),' out of ', num2str(Nframes)]);
    ii=ii+1;
end
disp('done')
% %%
% idxs=1:size(stack,3);
% for k=1:length(idxs)
%     imagesc(stack(:,:,idxs(k)));
%     colormap gray
%     title(files(k).name)
% %     waitforbuttonpress;
%     pause(0.005);
% end
