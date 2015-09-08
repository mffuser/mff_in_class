% script 1
x = [2 3; 4 2];
b = [1; 2];
bMatr = b(1:2, [1 1]);
c = x + bMatr;
c = x - bMatr;

matrProd = x*c;
matrProd/c

matrProd*inv(c)

%% elementwise multiplication

% initialize sale vector
nSold = [2 3 2 1];

prices = [2.3 2.4 5.4 1.2];

profits = nSold .* prices;

%% evaluate function

grid = 1:10;

squared = grid.^2;
squaredScaled = squared * 3;
scaled = 7*grid;
vals = squaredScaled + scaled

%% plotting

plot(grid, vals)

%%

grid = -4:0.02:4;
vals = exp(grid);

plot(grid, vals)

%%










