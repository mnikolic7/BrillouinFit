clear;
FOLDER='Oct23_m4_d5\';
FILENAMES=dir(FOLDER);
qq=1;
%
for nn=1:length(FILENAMES)
    if (~FILENAMES(nn).isdir)
%         [FOLDER,FILENAMES(nn).name]
        load([FOLDER,FILENAMES(nn).name]);
        Oct23_m4_d5(qq)=cell_data;
        qq=qq+1;
    end
end
%
save([FOLDER(1:end-1),'.mat'],FOLDER(1:end-1));
disp('done');