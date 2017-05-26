% Script file: S.O.M complex functions tested: 
% including: 
%        - REPULSION test.
%        - INTENSITY test. 
%        - BOTH (repulsion+intensity) tests.
% 
%% SINGLE FRAME/CLUMP
% Initialisation

tidy;
[dn,ds] = loadnames('macros', chooseplatform);
indx=48;

load(strcat(dn,ds(1:end-1),'_GT',filesep,'man00',num2str(indx),'.mat'));
[XREAL,xxa] = readParseInput(strcat(dn,ds,filesep,'man00',num2str(indx),'.tif'));

X = dataBin>0;

whichClump = 2;

[nuclei,clumps] = simpleClumpsSegmentation(double(X));
[boundies, clumpStr] = getOverlappingClumpsBoundaries(clumps, nuclei);
[angleMatrix,thisBoundy, numPoints] = computeAngleMatrix(boundies,whichClump);
[candies, candieshandles] = computeCandidatePoints(angleMatrix, ...
    thisBoundy, 150, 7, 'max');

thisClump=clumpStr.overlappingClumps==whichClump;
thisClumpNoNuclei = bitxor(thisClump,imdilate(bitand(thisClump, nuclei),ones(3)));

regs = regionprops(thisClump, 'BoundingBox');

% Deal with inputData here
X2=XREAL(:,:,2);

%% NORMAL INPUT
levs = multithresh(X2(clumpStr.labelledClumps>0));
dataIntensities = (XREAL(:,:,2)>levs).*thisClumpNoNuclei.*XREAL(:,:,2);
%A = ones(size(X2));

[inputData, sizeinput] = somGetInputData(dataIntensities.*A,XREAL);

%% BOUNDING BOX INPUT
A = zeros(size(X2));
%A(187:(205+20),265:(273+20)) = 1;
amin = min(candies,[],1);
amax = max(candies,[],1);
A(amin(1)-20:(amax(1)+20), amin(2)-20:(amax(2)+20)) = 1;

dataIntensities = X2.*A.*thisClumpNoNuclei;
dataIntensities(A>0) = dataIntensities(A>0)./max(dataIntensities(A>0));
[inputData, sizeinput] = somGetInputData(dataIntensities,XREAL);

%% WHOLE CLUMP INPUT (NO NUCLEI)
dataIntensities = X2.*thisClumpNoNuclei;
dataIntensities(thisClumpNoNuclei>0) = dataIntensities(thisClumpNoNuclei>0)./max(dataIntensities(thisClumpNoNuclei>0));
[inputData, sizeinput] = somGetInputData(dataIntensities,XREAL);

%% Complex network topology
net1size = 100;
net1type = 'ring';
%net1pos = thisBoundy(round(linspace(1,size(thisBoundy,1),100)),:);

net1pos = [candies(1,:); 
          thisBoundy(round(linspace(candyloc(1),candyloc(2),(net1size/2-1))),:)
          candies(2,:);
          thisBoundy(round(linspace(candyloc(2),numPoints,(net1size/4-1))),:)
          thisBoundy(round(linspace(1,candyloc(1),(net1size/4))),:)];

numJunctions = size(candies, 1);
net2size = numJunctions*10;
net2type = 'ring';
[~, srtcandies] = sort(candieshandles.intensityLocations);

if numJunctions > 2
    sortedcandies = candies(srtcandies,:);
else
    % only two junctions detected
    c1 = candies(srtcandies(1),:);
    c2 = candies(srtcandies(2),:);
    cm = mean(candies,1);
    v = c2-c1;
    v2 = v*[0 -1;1 0];
    corth1 = cm - 0.5.*v2;
    corth2 = cm + 0.25.*v2;
    
    sortedcandies = [c1;corth1;c2;corth2];
    numJunctions = size(sortedcandies, 1);
    net2size = numJunctions*10;
end

segsize = net2size/numJunctions - 1;
t = linspace(0.1,0.9,segsize)';

net2pos = [];
for ix=1:numJunctions
    SEG = repmat(sortedcandies(ix,:), segsize,1) + ...
        t.*repmat(sortedcandies(wrapN(ix+1,numJunctions),:) - ...
        sortedcandies(ix,:),segsize,1);
    net2pos = [net2pos;sortedcandies(ix,:); SEG];
end

[OG] = somComplexNetworks(net1pos, net2pos, sortedcandies, true);

figure(99)
somScatterGraphPlot(inputData, [], OG);

%% DISCONECTED CIRCLES
net1size = 100;
net1type = 'ring';
%net1pos = thisBoundy(round(linspace(1,size(thisBoundy,1),100)),:);
candyloc = candieshandles.intensityLocations;

net1pos = thisBoundy(round(linspace(candyloc(1),candyloc(2),net1size)),:);
speedouter1 = 0.25.*ones(size(net1pos,1),1);
net2pos = [thisBoundy(round(linspace(candyloc(2),numPoints,(net1size/2-1))),:)
          thisBoundy(round(linspace(1,candyloc(1),(net1size/2))),:)];
speedouter2 = 0.25.*ones(size(net1pos,1),1);

numJunctions = size(candies, 1);
net2size = numJunctions*10;
net2type = 'ring';
[~, srtcandies] = sort(candieshandles.intensityLocations);

if numJunctions > 2
    sortedcandies = candies(srtcandies,:);
else
    % only two junctions detected
    c1 = candies(srtcandies(1),:);
    c2 = candies(srtcandies(2),:);
    cm = mean(candies,1);
    v = c2-c1;
    v2 = v*[0 -1;1 0];
    corth1 = cm - 0.25.*v2;
    corth2 = cm + 0.25.*v2;
    
    sortedcandies = [c1;corth1;c2;corth2];
    numJunctions = size(sortedcandies, 1);
    net2size = numJunctions*10;
end

segsize = net2size/numJunctions - 1;
t = linspace(0.1,0.9,segsize)';

innernetpos = [];
for ix=1:numJunctions
    SEG = repmat(sortedcandies(ix,:), segsize,1) + ...
        t.*repmat(sortedcandies(wrapN(ix+1,numJunctions),:) - ...
        sortedcandies(ix,:),segsize,1);
    innernetpos = [innernetpos;sortedcandies(ix,:); SEG];
end

innernet1 = innernetpos(net2size/2+1:end,:);
speedinner1 = ones(size(innernet1,1),1);
innernet2 = innernetpos(1:net2size/2,:);
speedinner2 = ones(size(innernet2,1),1);

net1pos = [net1pos;innernet1];
net2pos = [net2pos;innernet2];

[OG] = somComplexNetworks(net1pos, net2pos, sortedcandies, false);
OG.Nodes.Speed = [speedouter1;speedinner1;speedouter2;speedinner2(1:end-1)];

%% THREE NETWORKS
[~, srtcandies] = sort(candieshandles.intensityLocations);

% only two junctions detected
c1 = candies(srtcandies(1),:);
c2 = candies(srtcandies(2),:);
cm = mean(candies(1:2,:),1);

% Get two grids in bounding boxes
netsize = [8 8];
nettype = 'supergrid';
amin = min(candies,[],1);
amax = max(candies,[],1);
bbox1 = [cm -8 -8];
bbox2 = [cm 8 8];

pos1 = somGetNetworkPositions(nettype, netsize, bbox1);
pos2 = somGetNetworkPositions(nettype, netsize, bbox2);
pos3 = thisBoundy(round(linspace(1,numPoints,100)),2:-1:1);

subnets.G1 = somBasicNetworks(nettype, netsize, pos1(:,2:-1:1));
subnets.G1.Nodes.ids = subnets.G1.Nodes.ids + 1000;
subnets.G2 = somBasicNetworks(nettype, netsize, pos2(:,2:-1:1));
subnets.G2.Nodes.ids = subnets.G2.Nodes.ids + 2000;
subnets.G3 = somBasicNetworks('ring',100,pos3);
subnets.G3.Nodes.ids = subnets.G3.Nodes.ids + 3000;
subnets.G3.Nodes.Speed = 0.1.*subnets.G3.Nodes.Speed;

netnames = fieldnames(subnets);

OG = graph;
ids = [];
x = [];
y = [];
Speed = [];

for kx=1:length(netnames)
    thisG = getfield(subnets, netnames{kx});
    auxedgetable = thisG.Edges;
    auxedgetable.EndNodes = auxedgetable.EndNodes + numnodes(OG);
    
    OG = addedge(OG, auxedgetable);
    %nodeindx = numnodes(OG)+(1:numnodes(thisG));
    
    ids = [ids;thisG.Nodes.ids];
    x = [x;thisG.Nodes.x];
    y = [y;thisG.Nodes.y];
    Speed = [Speed;thisG.Nodes.Speed];
end

OG.Nodes.ids = ids;
OG.Nodes.x = x;
OG.Nodes.y = y;
OG.Nodes.Speed = Speed;

%% parameters of the network
thisalpha = 0.25;
thisalphadecay = 'none';
thisneighboursize = 8;
thisneighbourdecay = 'exp';
thissteptype = 'intensity';

options.maxiter = 1000;
options.alphazero = thisalpha;
options.alphadtype = thisalphadecay;
options.N0 = thisneighboursize;
options.ndtype = thisneighbourdecay;
options.staticfrontier = true;
options.steptype = thissteptype;
options.debugvar = true;
options.gifolder = strcat('.',filesep);
options.gifname = 'labdemo.gif';
% options.gifname = strcat('teststep',upper(thissteptype), '-alpha',...
%     upper(thisalphadecay),'.gif');

figure(99)
%pause;
%%
tic;
[G, nethandles] = somTrainingPlus(inputData, OG, options);
nethandles.time = toc;

%% SOME PLOTS

inputstr = input(strcat('Name of string:',32),'s');
%inputstr = 'Test-image';

figure(10)
subplot(121)
somScatterGraphPlot(inputData,[],G,OG);
title(inputstr);
subplot(122)
plotGraphonImage(X2,OG);
plotGraphonImage([],G);
legend('Original Network', 'Final after 1000 iterations');
%axis([190 330 110 270])
