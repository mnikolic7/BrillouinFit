path='E:\current_data\20201019_m1234_d1\shape bf\m3_4\z_stack_single_cell_0.5um z step\';

files=dir([path,'*.tdms']);
Nframes=length(files);
I1=get_image_from_tdms([path,files(1).name]);
I1=flipud(rot90(I1));

%
Nrows=size(I1,1)
Ncolumns=size(I1,2)
stack=zeros(Nrows,Ncolumns,Nframes);
imagesc(I1)
%%
I=get_image_from_tdms([path,fname]);
I=uint16(I);
%%
imshow(I,[]);
%%
imwrite(I,'10A_60x_day5_2.tif');
%%
imshow(x10AT_60x_day5_7,[])
