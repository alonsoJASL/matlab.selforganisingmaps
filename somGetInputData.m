function [inputData, sizeinput] = somGetInputData(workingdata, rgbinput)
%           GET INPUT DATA FROM MATRIX or GENERIC
%
%

if nargin < 2
    rgbinput = [];
elseif size(rgbinput,3)==3
    red = rgbinput(:,:,1);
    green = rgbinput(:,:,2);
else
    red = zeros(size(rgbinput));
    green = rgbinput;
end

if ischar(workingdata)
    switch lower(workingdata)
        case  'square'
            inputData = rand(1000,2);
        case 'cross'
            b = rand(1000,2);
            b = [b(:,1) (b(:,2)/4)];
            b = [b(:,1) (b(:,2)+0.375)];
            c = [b(:,2) b(:,1)];
            a = [b;c];
            a = a(randperm(size(a,1)),:);
            
            inputData = a;
            clear a  b  c;
            
        case 'ball'
            r=0.5;
            a = rand(1000, 2);
            c = repmat([0.5 0.5], 1000, 1);
            
            b = a-c;
            b = sqrt(b(:,1).^2 + b(:,2).^2 );
            
            a = a(b < r, :);
            
            inputData = a;
            clear a c r;
        case 'triangle'
            a = rand(1000,2);
            b = inpolygon(a(:,1), a(:,2), [0;0.5;1],[0;1;0]);
            
            inputData = a(b,:);
            
            clear a b;
    end
else
    [xx, yy] = find(workingdata>0);
    inputData = [yy xx];
    if ~isempty(rgbinput)
        %red(workingdata>0) = red(workingdata>0)./max(red(workingdata>0));
        green(workingdata>0) = green(workingdata>0)./max(green(workingdata>0));
        
        inputData = [inputData red(workingdata>0) green(workingdata>0)];
    end
    
end

if nargout > 1
    sizeinput = size(inputData,1);
end