function [G, nethandles] = somTraining(inputData, initnetwork, options)
%               S.O.M. TRAINING FOR FULL DATA
%
%
% INITNETWORK: Either a structure with {type,size,pos}, or a Graph with the
% prebuilt Graph.
%
% OPTIONS: A structure containing the fields:
%          {alphazero, alphadtype, N0, ndtype, maxiter, debugvar}
%
%

if nargin<3
    maxiter = 1000;
    alphazero = 0.25;
    alphadtype = 'none';
    N0 = 3;
    ndtype= 'none';
    debugvar = false;
    gifolder = strcat('.',filesep);
    gifname = strcat('gifTest',num2str(randi(1000)),'.gif');
else
    [alphazero, alphadtype, N0, ndtype, maxiter, debugvar, ...
        gifolder, gifname] = getOptions(options);
end

options.alphazero = alphazero;
options.alphadtype = alphadtype;
options.N0 = N0;
options.ndtype = ndtype;
options.maxiter = maxiter;
options.debugvar = debugvar;

if isstruct(initnetwork)
    G = somBasicNetworks(initnetwork.type,initnetwork.size,initnetwork.pos);
else
    G = initnetwork;
end
OG = G;

xx = inputData(:,1);
yy = inputData(:,2);
sizeinput = length(xx);

if debugvar==true
    figure(99)
    somScatterGraphPlot(inputData, [], OG);
    title(strcat('Iter:', 32, num2str(0),'/', num2str(maxiter)));
    
    H = gcf;
    f = getframe(H);
    [im,map] = rgb2ind(f.cdata,256,'nodither');
    clf;
    
    im(1,1,1,maxiter) = 0;
    %pause;
end

% TYPES OF ALPHA DECAY: 'none', 'linear', 'exp1', 'exp2'
alphavect = alphazero.*alphaDecay(alphadtype,maxiter);
% TYPES OF NEIGHBOUR DECAY: 'none', 'linear', 'exp'
neighvect = neighbourDecay(ndtype, N0, maxiter);

for k=1:maxiter
    
    ridx = randi(sizeinput);
    
    thisinput = [xx(ridx) yy(ridx)];
    thisalpha = alphavect(k);
    thisnumneighbors = neighvect(k);
    
    [G, nethandles] = somSimpleStep(G, thisinput, thisalpha, thisnumneighbors);
    
    % debug
    if debugvar==true
        figure(99)
        clf
        somScatterGraphPlot(inputData, thisinput, G, nethandles.prevG);
        plot(OG,'--.b','XData',OG.Nodes.x,'YData',OG.Nodes.y);
        plot(nethandles.winnernode(1), nethandles.winnernode(2),'*g');
        title(strcat('Iter:', 32, num2str(k),'/', num2str(maxiter)));
        disp([k thisnumneighbors thisalpha]);
        
        f = getframe(gcf);
        im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
        
        %pause(0.01);
    end
end

if debugvar == true
    outputname = strcat(gifolder, gifname);
    imwrite(im,map,outputname,'DelayTime',0.25, 'LoopCount',inf);
end

if nargout > 1
    nethandles.originalG = OG;
    nethandles.options = options;
end

end

function [alphaval] = alphaDecay(decaytype, maxiter)
kvect = 1:maxiter;
switch lower(decaytype)
    case 'linear'
        alphaval = 1-kvect./maxiter;
    case 'exp1'
        alphaval = 0.2 + exp(-kvect./maxiter);
    case 'exp2'
        alphaval = 0.2 + exp(-kvect.^2./maxiter);
    otherwise
        alphaval = ones(1,maxiter);
end
end
function [neighval] = neighbourDecay(decaytype, N0, maxiter)
kvect = 1:maxiter;
switch lower(decaytype)
    case 'linear'
        % N1(x) = ceil(3 (1 - x / 1000))
        neighval = ceil(N0.*(1-kvect./maxiter));
    case 'exp'
        neighval = N0.*ceil(0.3 + 0.7.*exp(-kvect.^2./3000));
    otherwise
        neighval = N0.*ones(maxiter,1);
end
end
function [alphazero, alphadtype, N0, ndtype, maxiter, ...
    debugvar, gifolder, gifname] = getOptions(s)
% Get options for this function
%
alphazero = 0.25;
alphadtype = 'none';
N0 = 3;
ndtype = 'none';
maxiter = 1000;
debugvar = false;
gifolder = strcat('.',filesep);
gifname = strcat('gifTest',num2str(randi(1000)),'.gif');

fnames = fieldnames(s);
for ix=1:length(fnames)
    name = fnames{ix};
    switch name
        case 'alphazero'
            alphazero = getfield(s, name);
        case 'alphadtype'
            alphadtype = getfield(s, name);
        case 'N0'
            N0 = getfield(s, name);
        case 'ndtype'
            ndtype = getfield(s, name);
        case 'maxiter'
            maxiter = getfield(s, name);
        case 'debugvar'
            debugvar = getfield(s, name);
        case 'gifolder'
            gifolder = getfield(s, name);
        case 'gifname'
            gifname = getfield(s, name);
        otherwise
            str = strcat('[somTraining].ERROR. Incorrect option', 32,...
                'selected, field:',32, upper(name), 32,'is NOT defined.');
            disp(str);
    end
end
end





