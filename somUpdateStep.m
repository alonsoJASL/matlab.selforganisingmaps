function [newG, nethandles] = somUpdateStep(OG, inputxyrg, options)
%           S.O.M (SINGLE STEP) UPDATE NETWORK
%
% Computes one step of the S.O.M algorithm based on the current network
% (already initialised) and the input on which the network will adapt to.
%
% USAGE:
%        [newG, nethandles] = somUpdateStep(OG, inputXY, thisalpha, numneighbors)
%
% INPUT:
%                       OG := Original Graph with the network topology
%                            defined on its edges and the node positions
%                            under: [OG.Nodes.x OG.Nodes.y].
%                inputcoor := Input to the network. For initial
%                            implementation, only positions (x,y) are
%                            contemplated, however, this could include more
%                            features such as image intensities, i.e.
%                            (x,y,I(x,y)).
%                 options := Structure that contains
% OUTPUT:
%                     newG := Updated Network.
%               nethandles := Structure with previous network, index of
%                           winning node and the neighbours affected by the
%                           update.
%

if nargin < 3
    thisalpha = 0.25;
    numneighbors = 3;
    staticfrontier = false;
    steptype = 'none';
else
    [thisalpha, numneighbors, staticfrontier, steptype] = getOptions(options);
end

options.thisalpha = thisalpha;
options.numneighbors = numneighbors;
options.staticfrontier = staticfrontier;
options.steptype = steptype;

inputxy = inputxyrg(1:2);

if length(inputxyrg)>2
    repulsionAtt = inputxyrg(4) - inputxyrg(3);
    greenattenuation = inputxyrg(4);
else
    repulsionAtt = 1;
    greenattenuation = 1;
end

% make a winner node with an index (windx) on the network.
nodespositions = [OG.Nodes.x OG.Nodes.y];
testVect = nodespositions - repmat(inputxy,size(nodespositions,1),1);
testVec = testVect(:,1).^2 + testVect(:,2).^2;

[~,windx] = min(testVec);

dist2winnode = distances(OG,windx);
neighidx = find(dist2winnode <= numneighbors);

alphaval = thisalpha.*ones(length(neighidx),1);
neighbourpos = [OG.Nodes.x(neighidx) OG.Nodes.y(neighidx)];

% Take the step
%alphastep = alphaval.*(repmat(inputxy,size(neighbourpos,1),1) - neighbourpos);

winnernode = [OG.Nodes.x(windx) OG.Nodes.y(windx)];
%alphastep = alphaval.*repmat(inputxy - winnernode, size(neighbourpos,1),1);
alphastep = thisalpha.*repmat(inputxy - winnernode, size(neighbourpos,1),1);

if staticfrontier == true
    speedatt = OG.Nodes.Speed(neighidx);
    alphastep = alphastep.*speedatt;
end

switch lower(steptype)
    case 'none'
        updatedneighbours = neighbourpos + alphastep;
    case 'intensity'
        updatedneighbours = neighbourpos + (greenattenuation+0.2).^4.*alphastep;
    case 'repulsion'
        if repulsionAtt < 1
            updatedneighbours = neighbourpos - 0.25.*alphastep;
        else
            updatedneighbours = neighbourpos + alphastep;
        end
    case 'both'
        if repulsionAtt < 1
            updatedneighbours = neighbourpos + 0.25.*repulsionAtt.*alphastep;
        else
            updatedneighbours = neighbourpos + repulsionAtt.*alphastep;
        end
    otherwise
        updatedneighbours = neighbourpos + alphastep;
end

newG = OG;

newG.Nodes.x(neighidx) = updatedneighbours(:,1);
newG.Nodes.y(neighidx) = updatedneighbours(:,2);

if nargout > 1
    nethandles.stepOptions = options;
    nethandles.prevG = OG;
    nethandles.winnodeidx = windx;
    nethandles.winnernode = [OG.Nodes.x(windx) OG.Nodes.y(windx)];
    nethandles.neighborsaffected = neighidx;
end

end

function [thisalpha, numneighbors, staticfrontier, steptype] = getOptions(s)
% Get options for this function
%
thisalpha = 0.25;
numneighbors = 3;
staticfrontier = false;
steptype = 'none';

fnames = fieldnames(s);
for ix=1:length(fnames)
    name = fnames{ix};
    switch name
        case 'thisalpha'
            thisalpha = getfield(s, name);
        case 'numneighbors'
            numneighbors = getfield(s, name);
        case 'staticfrontier'
            staticfrontier = getfield(s, name);
        case 'steptype'
            steptype = getfield(s, name);
        otherwise
            str = strcat('[somTraining].ERROR. Incorrect option', 32,...
                'selected, field:',32, upper(name), 32,'is NOT defined.');
            disp(str);
    end
end
end
