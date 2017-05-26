% 
tidy;

whichFrame = 48;
whichClump = 2;


alphazero = 0.25; 
alphadtype = 'linear';
numneigh = 3;
neighdtype = 'linear';

[dn] = loadnames('macros', chooseplatform);

foldername = strcat(dn,'RESULTS', filesep, 'SOM',filesep, 'FRAME', ...
    num2str(whichFrame), filesep);
filename = strcat('clump', num2str(whichClump), '_alpha', upper(alphadtype),...
    num2str(alphazero*100), '_N-', upper(neighdtype),num2str(numneigh),'.mat');

disp('Loading...');
disp(strcat(foldername, filename));
load(strcat(foldername, filename));

somScatterGraphPlot([],[],G,nethandles.originalG);