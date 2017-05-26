function [XData, YData] = somGetPositionsFromBoundingBox(nettype, dataIn, ...
    netsizes, isRandom)
%             S.O.M NODE POSITIONS FROM BOUNDING BOX
%
%

if nargin < 2
    errmessage = strcat('[somGetPositionsFromBoundingBox].ERROR: ',...
        'Need to specify NETTYPE and give either a BOUNDING BOX or a', 32,...
        'BINARY IMAGE with a clump.');
    disp(errmessage);
    help somGetPositionsFromBoundingBox;
    
    XData = [];
    YData = [];
    return;
else
    if length(dataIn)>4
        % dataIn = a bounding box.
        bbox = dataIn;
        ox = bbox(1);
        oy = bbox(2);
        subwinx = fix(bbox(3)/netx);
        subwiny = fix(bbox(4)/nety);
        
        kx=1;
        xx = zeros(netx*nety,1);
        yy = zeros(netx*nety,1);
        
        for ix=1:netx
            for jx=1:nety
                xx(kx) = ox + (ix-1)*subwinx + fix(subwinx/2);%randi(subwinx);
                yy(kx) = oy + (jx-1)*subwiny + fix(subwiny/2);%randi(subwiny);
                kx=kx+1;
            end
        end
        
        XData = xx;
        YData = yy;
    else
        regs = regionprops(dataIn, 'BoundingBox');
        bbox = regs.BoundingBox;
        
        
    end
        