% Own SOM
% 
% Note that:
% 1. Input is defined as a (NN x 4) matrix with the columns:
%        [x y F(x,y)] = [x y f1(x,y) f2(x,y)]
% where:
%        f1(x,y) = intensity of image in red channel at position (x,y).
%        f2(x,y) = intensity of image in green channel at position (x,y).
% 
% 2. Nodes will be implemented in two ways, depending on how to handle the
% values of F(x,y):
%        2.1. NODES 1: Equations only work with position input, and
%        intensity values for F are used as attractors or repulsors of the
%        neighbouring nodes. 
%        This means that the input will only use values of intensity that
%        'help' it achieve the shape it wants: i.e values of intensity
%        higher than the second highest value for multithresh(greenIMAGE,2).
% 
%        2.2. NODES 2: Equations use both position and intensity values to
%        find the minimum and update the network. Input involves all the
%        positions and intensity values of the entire clump.
%

%% INITIALISATION
tidy;

[dn,ds] = loadnames('macros',chooseplatform);
indx=48;

load(strcat(dn,ds(1:end-1),'_GT',filesep,'man00',num2str(indx),'.mat'));
[XREAL,xxa] = readParseInput(strcat(dn,ds,filesep,'man00',num2str(indx),'.tif'));
Xedge = edge(XREAL(:,:,2),'canny',[],1);

GT = dataBin;

X = dataBin>0;
Xr = X(:,:,1);
Xg = X(:,:,2);
[M,N] = size(Xg);

%% BUILD OWN S.O.M Network 

whichClump = 2;
mainthresh = 150;
offsetVar = 7;
statsfname = 'max';

[nuclei,clumps] = simpleClumpsSegmentation(double(X));
[boundies, clumpStr] = getOverlappingClumpsBoundaries(clumps, nuclei);
[angleMatrix,thisBoundy, numPoints] = computeAngleMatrix(boundies,whichClump);
[candies, candieshandles] = computeCandidatePoints(angleMatrix, ...
    thisBoundy, mainthresh, offsetVar, statsfname);

thisClump=clumpStr.overlappingClumps==whichClump;

regs = regionprops(thisClump, 'Area', 'BoundingBox');

xx = randi(regs.BoundingBox(3),1000,1) + fix(regs.BoundingBox(1));
yy = randi(regs.BoundingBox(4),1000,1) +fix(regs.BoundingBox(2));

% Find points within the 
pinClump = inpolygon(xx,yy,thisBoundy(:,2), thisBoundy(:,1));
xin = xx(pinClump);
yin = yy(pinClump);

% Take the first nsom-2 points from the inside. Add candidates.
nsom = 40;
xnet = xin(1:nsom-2);
ynet = yin(1:nsom-2);

figure(1)
clf
imagesc(thisClump);
hold on
scatter(xx,yy,'+m');
scatter(xin,yin,'bo');
scatter(xnet,ynet,'dc');

netnodes = [xnet ynet];

%% NETWORKS FTW
xx = randi(regs.BoundingBox(3),40,1) + fix(regs.BoundingBox(1));
yy = randi(regs.BoundingBox(4),40,1) +fix(regs.BoundingBox(2));

G = somBasicNetworks('grid', [8 5]);

G.Nodes.x = xx;
G.Nodes.y = yy;
plotBoundariesAndPoints(XREAL, [], [yy xx],'m+')
plot(G,'XData', xx, 'YData',yy)


%% BUILD INPUTS FROM CLUMP, THRESHOLD AND INTENSITIES.

levs = multithresh(XREAL(:,:,2),2);

dataIntensities = bitor((XREAL(:,:,2)>levs(2)).*thisClump, ...
    bwperim(thisClump)).*XREAL(:,:,2);

%dataIntensities = (XREAL(:,:,2)>levs(2)).*thisClump.*XREAL(:,:,2);

[xx, yy] = find(dataIntensities>0);
inputxxyy = find(dataIntensities>0); 

%testKij = diag(XREAL(xx,yy,1));
testKij = cat(2,diag(XREAL(xx,yy,1)),diag(XREAL(xx,yy,2)));

dataIn = [xx yy testKij];

% input positions.
inputxy = dataIn(:,1:2);

figure(2)
imagesc(dataIntensities);

%% GET MINIMUM for NODES 1 (only positions) 

[pointsmindist, mindistval] = getBoundariesIntersection(netnodes,inputxy);
