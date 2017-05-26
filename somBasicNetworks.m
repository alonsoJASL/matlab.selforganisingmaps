function [gsom] = somBasicNetworks(nettype, netsize, nodepositions)
%                    S.O.M BASIC NETWORK
%
% INPUT:
%              gsom := Graph class object with network.
%                  - gsom.Edges := edges on the topology.
%                  - gsom.Nodes.XData := node positions in X
%                  - gsom.Nodes.YData := node positions in Y
% OUTPUT:
%           nettype := {'grid', 'ring', 'line'}
%           netsize := Size of network. 'grid' can receive one or
%                      two parameters in the vector. 'ring' and 'line'
%                      uses one value.
%     nodepositions := positions of the same size
%

if nargin == 0
    % default gris of 5x5
    nettype='grid';
    netsize=5;
    nodepositions=[];
elseif nargin < 2
    netsize = 5;
    nodepositions = [];
elseif nargin < 3
    nodepositions = [];
end

nettype = lower(nettype);
switch nettype
    case {'grid','supergrid'}
        if length(netsize) == 1
            AdjM = zeros(netsize^2);
            M = netsize;
            N = netsize;
        else
            AdjM = zeros(netsize(1)*netsize(2));
            M = netsize(1);
            N = netsize(2);
        end
        
        for ix=1:(M*N-1)
            AdjM(ix,ix + N) = 1;
            if mod(ix,N)~=0
                AdjM(ix,ix + 1) = 1;
            end
        end
        
        if strcmp(nettype,'supergrid')
            for ix=1:(M*N-1)
                if mod(ix,N) ~= 1
                    AdjM(ix,ix + (N-1)) = 1;
                end
                if mod(ix,N)~= 0
                    AdjM(ix,ix + (N+1)) = 1;
                end
            end
        end
        
        AdjM = AdjM(:,1:M*N) + AdjM(:,1:M*N)';
        
    case {'ring','circle'}
        AdjM = zeros(netsize);
        for jx = 1:netsize-1
            AdjM(jx,jx+1) = 1;
            if jx==1
                AdjM(jx,netsize) = 1;
            end
        end
        
        AdjM = AdjM + AdjM';
    case {'line'}
        AdjM = zeros(netsize);
        for jx = 1:netsize-1
            AdjM(jx,jx+1) = 1;
        end
        AdjM = AdjM + AdjM';
    otherwise
        errmessage = strcat('[',mfilename,'].ERROR:',...
            'Wrong NETTYPE.');
        disp(errmessage)
        help somBasicNetworks
        gsom = [];
        return;
end

gsom = graph(AdjM);
gsom.Nodes.Properties.Description = nettype;

if ~isempty(nodepositions)
    gsom.Nodes.ids = (1:length(nodepositions(:,1)))';
    gsom.Nodes.x = nodepositions(:,1);
    gsom.Nodes.y = nodepositions(:,2);
    gsom.Nodes.Speed = ones(size(gsom.Nodes,1),1);
end