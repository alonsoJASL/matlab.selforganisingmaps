function [newG, nethandles] = somSimpleStep(OG, inputxy, thisalpha, numneighbors)
%           S.O.M UPDATE NETWORK (1 STEP)
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
% OUTPUT:
%                     newG := Updated Network.
%               nethandles := Structure with previous network, index of
%                           winning node and the neighbours affected by the
%                           update.
%

if nargin < 3
    thisalpha = 0.25;
    numneighbors = 3;
elseif nargin < 4
    numneighbors = 3;
end

if isempty(thisalpha)==1
    thisalpha = 0.25;
end
if isempty(numneighbors)==1
    numneighbors = 3;
end

% make a winner node with an index (windx) on the network.
nodespositions = [OG.Nodes.x OG.Nodes.y];
testVect = nodespositions - repmat(inputxy,size(nodespositions,1),1);
testVec = testVect(:,1).^2 + testVect(:,2).^2;

[~,windx] = min(testVec); 

dist2winnode = distances(OG,windx);
neighidx = find(dist2winnode <= numneighbors); 

alphaval = thisalpha.*ones(length(neighidx),1);%exp(-dist2winnode(neighidx).^2/(2*numneighbors));
neighbourpos = [OG.Nodes.x(neighidx) OG.Nodes.y(neighidx)];

% Take the step
%alphastep = alphaval.*(repmat(inputxy,size(neighbourpos,1),1) - neighbourpos);

winnernode = [OG.Nodes.x(windx) OG.Nodes.y(windx)];
alphastep = alphaval.*repmat(inputxy - winnernode, size(neighbourpos,1),1);

updatedneighbours = neighbourpos + alphastep;

newG = OG;

newG.Nodes.x(neighidx) = updatedneighbours(:,1);
newG.Nodes.y(neighidx) = updatedneighbours(:,2);

if nargout > 1
    nethandles.prevG = OG;
    nethandles.winnodeidx = windx;
    nethandles.winnernode = [OG.Nodes.x(windx) OG.Nodes.y(windx)];
    nethandles.neighborsaffected = neighidx;
end

