%   Compute betti 1 curve of 2D root projection image to measure complexity
%   of root architecture
%
%   Description: Input : binary image with white foreground and black background 
%
%                Output : N-by-M matrix, where N is the number of images. M
%                is the number of sets
%
%   Author.....: Mao Li (maoli0923@gmail.com)
%   Date.......: Oct. 2016

file=dir('*.png'); % load the image, the format .png can be changed
N=length(file);
M=80; % a changeable paremeter which means the number of discreted sets for each betti 1 curve.
for k=1:N
    I=imread(file(k).name);
    [a b]=find(I>0);
    II=zeros(ceil(size(I)*1.2)); % line 5-8: the images were cropped without margin, we put a margin to each image, this could be changed
    for i=1:length(a)
        II(a(i)+200,b(i)+200)=255; 
    end
    for j=1:M
    A=ones(size(II));
    se = strel('disk',j-1);
    I2=imdilate(II,se); % dilate the image
    a=find(I2>0);
    A(a)=0; % switch black and white pixel 
    [labeled,H1(k,j)]=bwlabel(A,8); % compute the number of connected components of background (= number of holes of foreground+1)
    end
end
H1=H1-1; % number of holes of foreground= number of connected components of background -1.
save H1
