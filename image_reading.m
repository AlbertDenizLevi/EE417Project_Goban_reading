clear;
close;
imagefiles = dir('*.png');      
nfiles = length(imagefiles);    % Number of files found

for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   images{ii} = rgb2gray(currentimage);
end

%%
figure;
nfiles=5;
for ii=1:nfiles
    subplot(nfiles,1,ii);
    imshow(uint8(images{ii}));
    title(imagefiles(ii).name);
end
%%

