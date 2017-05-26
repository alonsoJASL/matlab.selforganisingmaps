function [G] = somComplexNetworks(nodepos1, nodepos2, sortedjunctions, joinedbool)
%

% Requires:
% sortedjunctions to come sorted in order of appearance inside thisBoundy
%

if nargin < 4
    joinedbool = false;
end

net1size = size(nodepos1,1);
net2size = size(nodepos2,1);
numJunctions = size(sortedjunctions, 1);
segmentsize = net2size/numJunctions;

G1 = somBasicNetworks('ring',net1size,nodepos1);
G2 = somBasicNetworks('ring',net2size,nodepos2);

G1.Nodes.Speed = 0.25.*G1.Nodes.Speed;
auxedgetable = G2.Edges;
auxedgetable.EndNodes = G2.Edges.EndNodes + net1size;

indx2erase = zeros(numJunctions,1);
if joinedbool == true
    for ix=1:numJunctions
        QQ = [G1.Nodes.x G1.Nodes.y] - repmat(sortedjunctions(ix,:),net1size,1);
        qq = QQ(:,1).^2 + QQ(:,2).^2;
        [minVal,minJix] = min(qq);
        
        if minVal< 30
            idx2change = net1size + (ix-1)*segmentsize + 1;
            indx2erase(ix) = idx2change;
            auxedgetable.EndNodes(auxedgetable.EndNodes==idx2change) = minJix;
        end
    end
    
    G1 = addedge(G1,auxedgetable);
    G1.Nodes.x(net1size+1:end) = G2.Nodes.x;
    G1.Nodes.y(net1size+1:end) = G2.Nodes.y;
    G1.Nodes.Speed(net1size+1:end) = G2.Nodes.Speed;
    
    indx2erase(indx2erase==0) = [];
    G1 = rmnode(G1,indx2erase);
    
else
    G1 = addedge(G1,auxedgetable);
    G1.Nodes.x(net1size+1:end) = G2.Nodes.x;
    G1.Nodes.y(net1size+1:end) = G2.Nodes.y;
    G1.Nodes.Speed(net1size+1:end) = G2.Nodes.Speed;
end

G = G1;
nx = G.Nodes.x;
G.Nodes.x = G.Nodes.y;
G.Nodes.y = nx;