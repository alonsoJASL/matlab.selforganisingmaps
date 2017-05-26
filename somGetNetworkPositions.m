function [pos, pos2] = somGetNetworkPositions(nettype, netsize, boundingbox, candies)
%        S.O.M GET NETWORK POSITION BASED ON NETWORK TYPE
%
%

if nargin < 4
    candies = [];
end

switch nettype
    case {'grid', 'supergrid'}
        netx = netsize(1);
        nety = netsize(2);
        
        ox = boundingbox(1);
        oy = boundingbox(2);
        subwinx = boundingbox(3)/netx;
        subwiny = boundingbox(4)/nety;
        
        % for random points in the grid:
        randomPoints = false;
        
        indxX = reshape(meshgrid(1:netx,1:nety),[netx*nety 1]);
        indxY = repmat((1:nety)',netx,1);
        
        if randomPoints == true
            xx = ox + (indxX-1).*subwinx + randi(subwinx,netx*nety,1);
            yy = oy + (indxY-1).*subwiny + randi(subwiny,netx*nety,1);
        else
            xx = ox + (indxX-1).*subwinx + subwinx/2;
            yy = oy + (indxY-1).*subwiny + subwiny/2;
        end
        
    case {'line','ring', 'circle'}
        if isempty(candies)
            ox = boundingbox(1) + boundingbox(3)/2;
            oy = boundingbox(2) + boundingbox(4)/2;
            
            circrad = 0.25*min(boundingbox([3 4]));
            
            numPoints = netsize;
            xx = cos(linspace(0,2*pi,numPoints))'.*circrad;
            yy = sin(linspace(0,2*pi,numPoints))'.*circrad;
            
            xx = xx+ox;
            yy = yy+oy;
            
        else
            midpoint = mean(candies,1);
            
            ox = midpoint(1);
            oy = midpoint(2);
            
            circrad = norm(candies(1,:)-midpoint);
            
            numPoints = netsize;
            xx = cos(linspace(0,2*pi,numPoints))'.*circrad;
            yy = sin(linspace(0,2*pi,numPoints))'.*circrad;
            
            xx = xx+ox;
            yy = yy+oy;
        end
 
    otherwise
        errmessage = strcat('[',mfilename,'].ERROR:',...
            'Wrong NETTYPE.');
        disp(errmessage)
        xx = [];
        yy = [];
        return;
        
end

pos = [xx yy];

if nargout == 2
    pos=xx;
    pos2=yy;
end