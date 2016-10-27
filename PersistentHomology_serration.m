%   Compute Euler characteristic (EC) curve (a persistent homology approach) to measure leaf serrations
%
%   Description: Input : binary image with black foreground and white background 
%                ( can be white foreground and black background, just delete line 26-28 and set h=I)
%
%                Output : N-by-M matrix, where N is the number of images, M
%                is the number of level sets (discrete a curve into some
%                sets). Each row has M numbers representing the EC curve to
%                measure leaf serrations.
%
%   Author.....: Mao Li (maoli0923@gmail.com)
%   Date.......: Oct. 2016
%   Remark.....: We use Alessandro Mannini and Auralius Manurung's code to
%                compute chain code and elliptical fourier transform
%                approximation

file=dir('*.jpg'); % load the image, the format .jpg can be changed
step=0.4/79;       % a changeable parameter of the interval size between two sets
N=length(file);    % number of images 
M=112;             % number of sets (changeable)
threshold=-0.3:step:-0.3+(M-1)*step; % thresholds of the sets
EC=zeros(N,M);     % initial the output matrix   

for k=1:length(file)  
    I=imread(file(k).name);  % read image
    h=zeros(size(I));        
    a=find(I==0);
    h(a)=1;                  % switch black and white pixel 
    [start_r,start_c] = find(h,1,'first'); % get start point
    V = bwtraceboundary(h,[start_r start_c],'W',8,Inf,'counterclockwise'); % trace the contour
    V=V-repmat(mean(V),length(V),1); % translate the center to origin 
    VV=zeros(size(V)); % initial VV
    VV(:,1)=V(:,2); VV(:,2)=V(:,1); % switch x y coordinates due to the difference between matrix index and image index 
    cc=chaincode(VV); % This is Alessandro Mannini's code to compute the chain code
    U=fourier_approx(transpose([cc.code]), 5, 1000, 0);% This is Auralius Manurung's code to get the EFT approximation, '5' is the number of harmonics, '1000' is the number of points on the approximated curve 
    U=U-repmat(mean(U),length(U),1); % translate the center to origin
    M=norm(V,'fro')/sqrt(length(V)); % centroid size of the leaf contour
    V=V/M;U=U/M;  % normalize both leaf contour and approximated curve
    DD=pdist2(V,U); % distance between leaf contour and approximated curve
    fun=zeros(length(V),1);
    in = inpolygon(V(:,1),V(:,2),U(:,1),U(:,2));
    for j=1:length(V)
        if in(j)
            fun(j)=-min(DD(j,:)); % inside the contour assign negative sign
        else
            fun(j)=min(DD(j,:));
        end
    end   % fun is the signed distance function
    E=[1:length(V);[2:length(V) 1]]'; % edge between two neighbor points
    fe=zeros(size(E,1),1); 
    fe=max(fun(E)')';  % assign value to edge connecting two points which is the maximum values of these two points's values
    for i=1:length(threshold)
        v=length(find(fun<=threshold(i))); %given a threshold, get the vertices whose values are smaller than it
        e=length(find(fe<=threshold(i)));  %given a threshold, get the faces whose values are smaller than it
        EC(k,i)=v-e; % Euler characteritic = number of vertics - number of faces
    end
end

save EC
