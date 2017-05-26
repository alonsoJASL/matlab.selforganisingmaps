% Script file: Networks and grids
%

%% INITIALISATION

tidy;

% load some image and some clumps to work.
[dn,ds] = loadnames('macros',chooseplatform);
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

regs = regionprops(thisClump, 'BoundingBox');

%% MAKING A RECTANGULAR GRID ON thisClump's BOUNDING BOX

% Build network
nettypes = {'grid', 'ring', 'line'};
netsizes = {[8 10], 8, 8};

whichnet = 1;

G = somBasicNetworks(nettypes{whichnet}, netsizes{whichnet});

netx = netsizes{whichnet}(1);
nety = netsizes{whichnet}(2);

ox = regs.BoundingBox(1);
oy = regs.BoundingBox(2);
subwinx = fix(regs.BoundingBox(3)/netx);
subwiny = fix(regs.BoundingBox(4)/nety);

% for random points in the grid:
randomPoints = false;

indxX = reshape(meshgrid(1:netx,1:nety),[netx*nety 1]);
indxY = repmat((1:nety)',netx,1);

if randomPoints == true
    xx = ox + (indxX-1).*subwinx + randi(subwinx,netx*nety,1);
    yy = oy + (indxY-1).*subwiny + randi(subwiny,netx*nety,1);
else
    xx = ox + (indxX-1).*subwinx + fix(subwinx/2);
    yy = oy + (indxY-1).*subwiny + fix(subwiny/2);
end

G.Nodes.x = xx;
G.Nodes.y = yy;

%% MAKING A RING/LINE NETWORK 

% get circle points 
% centre of the circle
ox = regs.BoundingBox(1) + fix(regs.BoundingBox(3)/2);
oy = regs.BoundingBox(2) + fix(regs.BoundingBox(4)/2);

circrad = 0.25*min(regs.BoundingBox([3 4]));

numPoints = 100;
xx = cos(linspace(0,2*pi,numPoints))'.*circrad;
yy = sin(linspace(0,2*pi,numPoints))'.*circrad;

xx = xx+ox;
yy = yy+oy;

% Build network
nettypes = {'grid', 'ring', 'line'};
netsizes = {[8 5], numPoints, 10};

whichnet = 2;

G = somBasicNetworks(nettypes{whichnet}, netsizes{whichnet});

G.Nodes.x = xx;
G.Nodes.y = yy;

%% MAKING TOO SELF-CONTAINING RINGS 

% Get position points and draw extended arcs to both. 
candyloc = candieshandles.intensityLocations;

xini = candies(1,:);
xend = candies(2,:);

cl = linspace(0.1,0.9,10);
Xextini = zeros(10,2);
Xextend = zeros(10,2);

for ix=1:10
    Xextini(ix,:) = cl(ix).*xini + (1-cl(ix)).*xend + randi([2 3],1,2);
    Xextend(ix,:) = cl(ix).*xini + (1-cl(ix)).*xend - randi([2 3],1,2);
end

segment1 = thisBoundy(candyloc(1):candyloc(2),:);
segment1 = [segment1;Xextini];

% make rings individually.
G1 = somBasicNetworks('ring', 8);
G2 = somBasicNetworks('ring', 10);

% changing node IDs so they don't overlap.
auxedgetable = G2.Edges;
auxedgetable.EndNodes = G2.Edges.EndNodes + numnodes(G1);

% find 
find(auxedgetable.EndNodes==9)
auxedgetable.EndNodes([1 2]) = 4;
find(auxedgetable.EndNodes==16)
auxedgetable.EndNodes([9 18]) = 1;
G1 = addedge(G1,auxedgetable);
G1 = rmnode(G1,[9 16]);

% add edges so rings are connected.

%% MAKING TOO INTERSECTING RINGS WITH CELL POINTS

% Get position points and draw extended arcs to both. 
candyloc = candieshandles.intensityLocations;

xini = candies(1,:);
xend = candies(2,:);

cl = linspace(0.1,0.9,10);
Xextini = zeros(10,2);
Xextend = zeros(10,2);

for ix=1:10
    Xextini(ix,:) = cl(ix).*xini + (1-cl(ix)).*xend + randi([2 3],1,2);
    Xextend(ix,:) = cl(ix).*xini + (1-cl(ix)).*xend - randi([2 3],1,2);
end

segment1 = thisBoundy(candyloc(1):candyloc(2),:);
segment1 = [segment1;Xextini];

% make rings individually.
G1 = somBasicNetworks('ring', 8);
G2 = somBasicNetworks('ring', 10);

% changing node IDs so they don't overlap.
auxedgetable = G2.Edges;
auxedgetable.EndNodes = G2.Edges.EndNodes + numnodes(G1);

% find 
find(auxedgetable.EndNodes==9)
auxedgetable.EndNodes([1 2]) = 4;
find(auxedgetable.EndNodes==16)
auxedgetable.EndNodes([9 18]) = 1;
G1 = addedge(G1,auxedgetable);
G1 = rmnode(G1,[9 16]);

% add edges so rings are connected.


%% UPDATE WINNER NODE

figure(1)
clf;
plotBoundariesAndPoints(XREAL,thisBoundy); hold on;
plot(G,'XData', G.Nodes.x, 'YData',G.Nodes.y);

inputXY = [283 184];
plot(inputXY(1), inputXY(2), 'm+');

% make a proper winnode
nodesIX = [G.Nodes.x G.Nodes.y]; % this iteration's nodes
testVect = nodesIX - repmat(inputXY,size(nodesIX,1),1);
testVec = sqrt(testVect(:,1).^2 + testVect(:,2).^2);

alphaIX = 0.5;

% calculate an artificial alphastep
[~,windx] = min(testVec); % winner node index m_c(t) in paper.
winnode = [G.Nodes.x(windx(1)) G.Nodes.y(windx(1))];
neighdist = 2; % Neighborhood distance.
winneighidx = find(distances(G,windx(1))<=neighdist); % N_c(t) in paper

% get neighborhood, update it and the Neighborhood in the graph
neighXY = [G.Nodes.x(winneighidx) G.Nodes.y(winneighidx)];
alphastep = alphaIX.*(repmat(inputXY,size(neighXY,1),1) - neighXY);

newNeighXY = neighXY + alphastep;

G.Nodes.x(winneighidx) = newNeighXY(:,1);
G.Nodes.y(winneighidx) = newNeighXY(:,2);

pause;
plot(G,'XData', G.Nodes.x, 'YData',G.Nodes.y);


%% VISUALISATION TOOL FOR INPUT DATA

% alpha decay: alphaIX = 

for indx=1:10
    figure(100)
    clf;
    plotBoundariesAndPoints(XREAL,thisBoundy); hold on;
    axis([xmmin xmmax ymmin ymmax]);
    plot(G,'XData', G.Nodes.x, 'YData',G.Nodes.y);
    
    [qx, qy] = ginput(1);
    inputXY = [qx, qy];
    plot(inputXY(1), inputXY(2), 'm+');
    
    % make a proper winnode
    nodesIX = [G.Nodes.x G.Nodes.y]; % this iteration's nodes
    testVect = nodesIX - repmat(inputXY,size(nodesIX,1),1);
    testVec = sqrt(testVect(:,1).^2 + testVect(:,2).^2);
    
    alphaIX = linspace(0.8,0.01,10);
    
    % calculate an artificial alphastep
    [~,windx] = min(testVec); % winner node index m_c(t) in paper.
    winnode = [G.Nodes.x(windx(1)) G.Nodes.y(windx(1))];
    neighdist = fix(linspace(mean(DD)/2,2,10));%[6 5 5 4 4 3 3 2 2 1]; % Neighborhood distance.
    winneighidx = find(distances(G,windx(1))<=neighdist(indx)); % N_c(t) in paper
    
    % get neighborhood, update it and the Neighborhood in the graph
    neighXY = [G.Nodes.x(winneighidx) G.Nodes.y(winneighidx)];
    alphastep = alphaIX(indx).*(repmat(inputXY,size(neighXY,1),1) - neighXY);
    
    newNeighXY = neighXY + alphastep;
    
    G.Nodes.x(winneighidx) = newNeighXY(:,1);
    G.Nodes.y(winneighidx) = newNeighXY(:,2);
    
    plot(G,'XData', G.Nodes.x, 'YData',G.Nodes.y);
    pause(0.5);
end

