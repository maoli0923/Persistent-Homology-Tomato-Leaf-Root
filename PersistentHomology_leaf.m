%   Compute persisetent barcodes for 16 rings for each leaf and pairwise bottleneck distance (a
%   persistent homology approach) to measure leaf shape
%
%   Description: Input : binary image with black foreground and white background 
%                ( can be white foreground and black background, just delete line 54-56 and set h=I)
%
%                Output : N-by-N matrix, where N is the number of images.
%                (i,j)-element is the distance between leaf i and leaf j.
%                Use multidimensional scaling technique, this matrix could
%                be projected into lower dimensional Euclidean space
%                similar to PCA.
%
%   Author.....: Mao Li (maoli0923@gmail.com)
%   Date.......: Oct. 2016
%   Remark.....: We use Javaplex to compute persistent barcodes and bottleneck distance

%Step 1, disrete the subset of 2D plane which is large enough to include leaf into simplical complex (some triangles).
%Those could be changed, but we use them for this study.
bound1=-2;bound2=2; % set boundary value
extra=0; % margin parameter
l1=100;  % number of sets in x axis
l2=100;  % number of sets in y axis
hgrid=bound1-extra*(bound2-bound1)/l1:(bound2-bound1)/l1:bound2+extra*(bound2-bound1)/l1; %discrete x axis
vgrid=hgrid; % discrete y axis which is same as x axis here
[ysample,xsample]=meshgrid(vgrid,hgrid); % grid in 2D plane
U(:,1)=reshape(xsample,size(xsample,1)*size(xsample,2),1);U(:,2)=reshape(ysample,size(ysample,1)*size(ysample,2),1); %re-format, U is the vertices on 2D plane meshes
n1=length(hgrid);n2=length(vgrid); 
E=[];
for i=1:n2-1
    for j=(i-1)*n1+1:i*n1-1
        E=[E;j j+1];
        E=[E;j j+n1];
        E=[E;j j+n1+1];
    end
    j=i*n1;
    E=[E;j j+n1];
end
i=n2;
for j=(i-1)*n1+1:i*n1-1
    E=[E;j j+1];
end  % E is the edges on 2D plane meshes

% Step 2, Compute persistent homology
% (Need to install javaplex package in Matlab)

% Those are parameters that could be changed. 
sigma0=0.085; %bandwith of Gaussian density estimator 
sigma=0.1;    %parameter to control the size of annulus
kN=16;        %number of rings
    
file=dir('*.jpg');
for k=1:length(file)
    I=imread(file(k).name); % load the image, the format .jpg can be changed
    h=zeros(size(I));
    a=find(I==0);
    h(a)=1;                 % switch black and white pixel 
    [start_r,start_c] = find(h,1,'first'); % get start point
    V = bwtraceboundary(h,[start_r start_c],'W',8,Inf,'counterclockwise'); % trace the contour
    V=V-repmat(mean(V),length(V),1); % translate the center to origin 
    V=V/norm(V,'fro')*sqrt(length(V)); %normalize leaf contour by centroid size
    F = ksdensity2d([V(:,1),V(:,2)],hgrid,vgrid,[sigma0,sigma0]); % Gaussian density estimator, some other person's code
    for t=1:kN 
        G=exp(-abs(sum((U-repmat([0 0],size(U,1),1)).^2,2).^0.5-t*sigma).^2/(2*sigma^2)); % annulus kernel
        Z=reshape(reshape(F,size(G)).*G,size(F)); 
        ZZ=-reshape(Z,n1*n2,1); % negative function so that swich superlevel set to sublevel set
        % The following is to compute persistent barocdes for each ring
        stream = api.Plex4.createExplicitSimplexStream(0.001); 
        for i=1:length(U)
            stream.addVertex(i,ZZ(i));
        end
        for i=1:length(E)
            stream.addElement([E(i,1) E(i,2)],max(ZZ(E(i,1)),ZZ(E(i,2))));
        end
        stream.finalizeStream();
        persistence = api.Plex4.getModularSimplicialAlgorithm(1, 2);
        intervals = persistence.computeAnnotatedIntervals(stream);
        Interval0{k}{t}=intervals.getIntervalsAtDimension(0); 
    end
end

%Step 3, compute pairwise bottleneck distance
BD=zeros(length(file),length(file),kN); 
for t=1:kN
    for i=1:size(BD,1)-1
        for j=i+1:size(BD,2)
            BD(i,j,t) = edu.stanford.math.plex4.bottleneck.BottleneckDistance.computeBottleneckDistance(Interval0{i}{t},Interval0{j}{t});
        end
    end
end
BD_leaf=zeros(length(file),length(file));
for t=1:kN
    BD_leaf=BD_leaf+BD(:,:,t).^2;
end
BD_leaf=sqrt(BD_leaf);
BD_leaf=BD_leaf+BD_leaf'; % BD_leaf is the final pairwise distance
save BD_leaf
