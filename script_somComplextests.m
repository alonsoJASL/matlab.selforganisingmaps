% script file: Complex S.O.M for overlapping segmentation 
% 

tidy;
whichCellTest = 'COMPLEXSOM';
[dn,ds] = loadnames('macros', chooseplatform);
matfolder = strcat('RESULTS',filesep,whichCellTest,filesep);

indx=50;
for indx=[45 46 48 49 50 51 52 53 54]
framefolder = strcat('FRAME',num2str(indx),filesep);

% Create output folders.
thisfoldername = strcat(dn,matfolder,framefolder);

if ~isdir(thisfoldername)
    mkdir(thisfoldername);
end

load(strcat(dn,ds(1:end-1),'_GT',filesep,'man00',num2str(indx),'.mat'));
[XREAL,xxa] = readParseInput(strcat(dn,ds,filesep,'man00',num2str(indx),'.tif'));
X2 = XREAL(:,:,2);

X = dataBin>0;
[nuclei,clumps] = simpleClumpsSegmentation(double(X));
[boundies, clumpStr] = getOverlappingClumpsBoundaries(clumps, nuclei);

imagesc(clumpStr.overlappingClumps)
[yq, xq] = ginput(1);

whichClump = clumpStr.overlappingClumps(fix(xq),fix(yq));
close all

[angleMatrix,thisBoundy, numPoints] = computeAngleMatrix(boundies,whichClump);
[candies, candieshandles] = computeCandidatePoints(angleMatrix, ...
    thisBoundy, 150, 7, 'max');

thisClump=clumpStr.overlappingClumps==whichClump;
thisClumpNoNuclei = bitxor(thisClump,imdilate(bitand(thisClump, nuclei),ones(3)));

% INPUT
A = zeros(size(X2));
%A(187:(205+20),265:(273+20)) = 1;
amin = min(candies,[],1);
amax = max(candies,[],1);
A(amin(1)-20:(amax(1)+20), amin(2)-20:(amax(2)+20)) = 1;

levs = multithresh(X2(clumpStr.labelledClumps>0));
dataIntensities = (XREAL(:,:,2)>levs).*thisClumpNoNuclei.*XREAL(:,:,2);

[inputData, sizeinput] = somGetInputData(dataIntensities.*A,XREAL);

% THREE NETWORKS
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

% parameters of the network
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
options.debugvar = false;
options.gifolder = strcat('.',filesep);
options.gifname = [];

tic;
[G, nethandles] = somTrainingPlus(inputData, OG, options);
nethandles.time = toc;

% Now save everything
thismatname = strcat('clump', num2str(whichClump),'_',lower(whichCellTest),...
    '-alpha0',num2str(thisalpha*1000),upper(thisalphadecay),...
    '-N',num2str(thisneighboursize),'-steptype',upper(thissteptype),'.mat');
%
save(strcat(thisfoldername, thismatname),...
    'thisClump','candies', 'candieshandles', 'OG', 'G', 'nethandles', ...
    'inputData', 'XREAL');
if strcmp(chooseplatform,'linux')
    onedrive = strcat('/home/jsolisl/OneDrive/MACROPHAGES/RESULTS/',...
        whichCellTest, filesep);
    framedrive = strcat(onedrive,framefolder);
    if ~isdir(framedrive)
        mkdir(framedrive);
    end
    
    save(strcat(framedrive,thismatname),...
        'thisClump','candies', 'candieshandles', 'OG', 'G', 'nethandles', ...
    'inputData', 'XREAL');
end
end

%%
inputstr = strcat('Test-image-');

figure(10)
subplot(121)
somScatterGraphPlot(inputData,[],G,OG);
title(inputstr);
subplot(122)
plotGraphonImage(X2,OG);
plotGraphonImage([],G);
legend('Original Network', 'Final after 1000 iterations');