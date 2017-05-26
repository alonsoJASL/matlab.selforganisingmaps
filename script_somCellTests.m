% script file: COMPREHENSIVE CELL TESTS
%
% Cell tests for many parameters scenarios such as:
%
% (image) indx = [45:54];
% alphazero = [0.125, 0.25];
%             - with decay types: {'none','linear','exp1', 'exp2'}
% N0 = [3 4];
%             - with decay types: {'none', 'linear', 'exp'}
%
% Organisation of the RESULTS of these experiments:
%   PATH_TO_/MACROPHAGES/RESULTS/SOM/FRAME<indx>/
%
% GIFs of the results will be saved on:
%  PATH_TO_/MACROPHAGES/GIF/SOM/FRAME<indx>/
%
% With file names (whether .mat or .gif):
%  'clump<whichClump>_alpha<DECAY><100alphazero>_N<DECAY><N0>.mat'
%

%% SINGLE FRAME/CLUMP
% Initialisation

tidy;
[dn,ds] = loadnames('macros', chooseplatform);
matfolder = strcat('RESULTS',filesep,'SOM',filesep);
giffolder = strcat('GIF',filesep,'SOM',filesep);

indx=48;
framefolder = strcat('FRAME',num2str(indx),filesep);

% Create output folders.
thisfoldername = strcat(dn,matfolder,framefolder);
thisgifolder = strcat(dn,giffolder);
if ~isdir(thisfoldername)
    mkdir(thisfoldername);
end
if ~isdir(thisgifolder)
    mkdir(thisgifolder);
end

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

% Deal with inputData here
levs = multithresh(XREAL(:,:,2),2);
dataIntensities = (XREAL(:,:,2)>levs(2)).*thisClumpNoNuclei.*XREAL(:,:,2);
[inputData, sizeinput] = somGetInputData(dataIntensities);

% Complex network topology
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

% parameters of the network
thisalpha = 0.25;
thisalphadecay = 'none';
thisneighboursize = 3;
thisneighbourdecay = 'none';

thisoutputname = strcat('clump',num2str(whichClump),'_alpha',...
    upper(thisalphadecay),num2str(thisalpha*100),'_N-',...
    upper(thisneighbourdecay), num2str(thisneighboursize));
thismatname = strcat(thisoutputname,'.mat');

options.maxiter = 1000;
options.alphazero = thisalpha;
options.alphadtype = thisalphadecay;
options.N0 = thisneighboursize;
options.ndtype = thisneighbourdecay;
options.debugvar = false;
options.gifolder = thisgifolder;
options.gifname = strcat(thisoutputname,'.gif');

tic;
[G, nethandles] = somTraining(inputData, OG, options);
nethandles.time = toc;

save(strcat(thisfoldername,thismatname),'G','nethandles');


%% ALL FRAMES/CLUMPS
% Initialisation

tidy;
[dn,ds] = loadnames('macros', chooseplatform);
matfolder = strcat('RESULTS',filesep,'SOM',filesep);
giffolder = strcat('GIF',filesep,'SOM',filesep);

allalphazeros = 0.125:0.025:0.25;
allalphadecay = {'none','linear','exp1', 'exp2'};
allneighboursizes = 2:4;
allneighbourdecay = {'none', 'linear', 'exp'};


for indx=45:54
    framefolder = strcat('FRAME',num2str(indx),filesep);
    
    % Create output folders.
    thisfoldername = strcat(dn,matfolder,framefolder);
    thisgifolder = strcat(dn,giffolder);
    if ~isdir(thisfoldername)
        mkdir(thisfoldername);
    end
    if ~isdir(thisgifolder)
        mkdir(thisgifolder);
    end
    
    load(strcat(dn,ds(1:end-1),'_GT',filesep,'man00',num2str(indx),'.mat'));
    [XREAL,xxa] = readParseInput(strcat(dn,ds,filesep,'man00',num2str(indx),'.tif'));
    
    X = dataBin>0;
    
    [nuclei,clumps] = simpleClumpsSegmentation(double(X));
    [boundies, clumpStr] = getOverlappingClumpsBoundaries(clumps, nuclei);
    
    nOverlapClumps = unique(clumpStr.overlappingClumps);
    for whichClump=1:max(nOverlapClumps)
        [angleMatrix,thisBoundy, numPoints] = computeAngleMatrix(boundies,whichClump);
        [candies, candieshandles] = computeCandidatePoints(angleMatrix, ...
            thisBoundy, 150, 7, 'max');
        
        thisClump=clumpStr.overlappingClumps==whichClump;
        thisClumpNoNuclei = thisClump-bitand(thisClump, nuclei);
        
        regs = regionprops(thisClump, 'BoundingBox');
        
        % Deal with inputData here
        levs = multithresh(XREAL(:,:,2),2);
        dataIntensities = (XREAL(:,:,2)>levs(2)).*thisClumpNoNuclei.*XREAL(:,:,2);
        [inputData, sizeinput] = somGetInputData(dataIntensities);
        
        % Complex network topology
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
        
        % parameters of the network
        for ax=1:length(allalphazeros)
            for dax=1:length(allalphadecay)
                for nx=1:length(allneighboursizes)
                    for dnx=1:length(allneighbourdecay)
                        
                        thisalpha = allalphazeros(ax);
                        thisalphadecay = allalphadecay{dax};
                        thisneighboursize = allneighboursizes(nx);
                        thisneighbourdecay = allneighbourdecay{dnx};
                        
                        thisoutputname = strcat('clump',num2str(whichClump),'_alpha',...
                            upper(thisalphadecay),num2str(thisalpha*100),'_N-',...
                            upper(thisneighbourdecay), num2str(thisneighboursize));
                        thismatname = strcat(thisoutputname,'.mat');
                        
                        options.maxiter = 1000;
                        options.alphazero = thisalpha;
                        options.alphadtype = thisalphadecay;
                        options.N0 = thisneighboursize;
                        options.ndtype = thisneighbourdecay;
                        options.debugvar = true;
                        options.gifolder = thisgifolder;
                        options.gifname = strcat(thisoutputname,'.gif');
                        
                        tic;
                        [G, nethandles] = somTraining(inputData, OG, options);
                        nethandles.time = toc;
                        
                        save(strcat(thisfoldername,thismatname),'G','nethandles');
                    end
                end
            end
        end
    end
end
