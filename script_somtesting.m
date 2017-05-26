% Script file: The definitive script to test my SOMs
%

%% INITIALISATION
tidy;

switch 'generic'
    case 'generic'
        inputtype = 'cross'; % square, cross, ball, triangle, donut
        [inputData, sizeinput] = somGetInputData(inputtype);
        regs.BoundingBox = [2 2 4 4];
        
    case {'ellipses', 'ellipse'}
        fname = '40_dist_30_ecc_90.mat'; % FOUR CANDIDATES
        %fname = '40_dist_80_ecc_90.mat'; % TWO CANDIDATES
        [dn,ds] = loadnames('macrosynth',chooseplatform);
        pathname = strcat(dn,ds(1:end-1),'_OriginalMax',filesep,fname);     
        
        
        load(pathname);
        dataIntensities = edge(X,'Canny');
        regs = regionprops(X>0,'BoundingBox');
        [inputData, sizeinput] = somGetInputData(dataIntensities);
        
        XREAL = X;
        
    case {'cells', 'cell'}
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
        thisClumpNoNuclei = thisClump-bitand(thisClump, nuclei);
        
        regs = regionprops(thisClump, 'BoundingBox');
        
        levs = multithresh(XREAL(:,:,2),2);
        dataIntensities = (XREAL(:,:,2)>levs(2)).*thisClumpNoNuclei.*XREAL(:,:,2);
        
        % add boundary perimeter (to keep things 'tidy')
        %dataIntensities = bitor(dataIntensities>0,bwperim(thisClump));
        
        [inputData, sizeinput] = somGetInputData(dataIntensities);
end

%% (CLASSIC) NETWORK TOPOLOGIES

% choose initial size of network.
%boundingbox = regs.BoundingBox;
boundingbox = [1 1 0.5 0.5];

netstruct.type = 'supergrid';
if strcmp(netstruct.type,'supergrid')
    netstruct.size = [10 10];
else
    netstruct.size = 200;
end

netstruct.pos = somGetNetworkPositions(netstruct.type,netstruct.size,...
    boundingbox);
OG = somBasicNetworks(netstruct.type,netstruct.size,netstruct.pos);

figure(1)
clf;
title('Data and initial network')
somScatterGraphPlot(inputData, [], OG);

%% TRAINING initial PARAMETERS

k=1;
nodedist = distances(OG);
N0 = max(nodedist(~isinf(nodedist)));
maxiter = 1000;

clc;
%% SINGLE STEP (FROM CHAOS)
krand = randi(sizeinput);
thisinput = [inputData(krand,1) inputData(krand,2)];

thisalpha = 0.25;
if k < maxiter/2
    numneighbors = 9 + N0*exp(-k/N0);
else
    numneighbors = 4 + N0*exp(-k/N0);
end
[OG, nethandles] = somSimpleStep(OG, thisinput, thisalpha, numneighbors);
disp([k numneighbors thisalpha]);

%if mod(k,100)==0
figure(3)
clf;
title(strcat('Iteration = ',32,num2str(k)));
somScatterGraphPlot(inputData, thisinput, OG, nethandles.prevG);
%end

k=k+1;

%% (COMPLEX) NETWORK TOPOLOGY

% ONLY with ellipses and cells
net1size = 100;
net1type = 'ring';
net1pos = thisBoundy(round(linspace(1,size(thisBoundy,1),100)),:);

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
    corth1 = cm - 0.125.*v2;
    corth2 = cm + 0.125.*v2;
    
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
OG = OG;

%[inputData, sizeinput] = somGetInputData(edge(dataIntensities));

figure(2)
clf;
title('Data and initial network (complex topology)')
somScatterGraphPlot(inputData, [], OG);


%% TRAINING initial PARAMETERS

k=1;
nodedist = distances(OG);
N0 = round(max(nodedist(~isinf(nodedist)))/2);
maxiter = 1000;

clc;
%% SINGLE STEP (TWO CIRCLES)
krand = randi(sizeinput);
thisinput = [inputData(krand,1) inputData(krand,2)];


thisalpha = 0.25;
numneighbors = 3;

stepOptions.thisalpha = thisalpha;
stepOptions.numneighbors = numneighbors;
stepOptions.steptype = 'none';
%%
[OG, nethandles] = somUpdateStep(OG,thisinput)
%[G, nethandles] = somSimpleStep(G, thisinput, thisalpha, numneighbors);
disp([k numneighbors thisalpha]);

%if mod(k,100)==0
    figure(3)
    clf;
    title(strcat('Iteration = ',32,num2str(k)));
    somScatterGraphPlot(inputData, thisinput, OG, nethandles.prevG);
    plot(OG,'--.b','XData',OG.Nodes.x,'YData',OG.Nodes.y);
    plot(nethandles.winnernode(1), nethandles.winnernode(2),'*g');
    title(strcat('ITER=',num2str(k),32,'MORE POINTS IN NETWORK'));
    axis([185 330 115 270]);
%end

k=k+1;

%% FULL TRAINING 

clear options;

options.maxiter = 100;
options.alphazero = 0.125;
options.alphadtype = 'linear';
options.N0 = 5;
options.ndtype = 'exp';
options.debugvar = true;

[Gfull, fullnethandles] = somTraining(inputData, OG, options);

%%
figure(4)
clf;
subplot(121)
somScatterGraphPlot(inputData, [],Gfull,fullnethandles.originalG);
plot(OG,'--.b','XData',OG.Nodes.x,'YData',OG.Nodes.y);

grid on;
subplot(122)
plotBoundariesAndPoints(XREAL,thisBoundy,candies);
plotGraphonImage([],Gfull,'+-m');
legend('','Clump boundary','junctions','S.O.M');
axis off;