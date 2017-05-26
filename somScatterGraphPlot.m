function somScatterGraphPlot(inputData, somepoint, G, OG)
% 

if nargin < 4
    OG =[];
end

if ~isempty(inputData)
    scatter(inputData(:,1), inputData(:,2));
end
hold on;
if ~isempty(OG)
    plot(OG,'+--r','XData', OG.Nodes.x, 'YData',OG.Nodes.y, ...
        'LineWidth',3);
end
plot(G,'dk','XData', G.Nodes.x, 'YData',G.Nodes.y,...
        'LineWidth',3);
if ~isempty(somepoint)
    plot(somepoint(1), somepoint(2),'dy', 'markersize',4);
end
grid on
axis ij
if nargin == 4
    legend('Data','Previous Net', 'Current Net', 'Input');
end
