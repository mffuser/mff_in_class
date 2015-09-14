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

grid = -4:0.00002:4;

mu=0;
sigma=1;

tic;
y = zeros(length(grid), 1);

for ii=1:length(grid)
    y(ii, 1) = (2*pi*sigma)^(-0.5)*exp(-0.5*((grid(ii)-mu)^2/sigma));
end

t = toc;

%%

y = zeros(3, 1);

for ii=1:3
    y(ii) = ii
end

%%

grid = -4:0.02:4;
y = (2*pi*sigma)^(-0.5) * exp(-0.5*(grid-mu).^2/sigma);

%%

a = [ 1  2  1 15;
      0  8  5  9];
  
most_change(a)

%%

rank(a(:,1))

%%

x = [1 2 3; 4 1 5; 0 4 4];
y_correct = [2 2 1; 3 1 3; 1 3 2];
isequal(ranks(x),y_correct)


%%

x = rand(3)

transp = x'

%%

x(:)'

%%

mu = [4; 0];
sigma = [1.5; 1];
rho = 0.5;

grid = -4:0.1:8;

tic

% preallocation
nGrid = numel(grid);
z = zeros(nGrid, nGrid);
xGrid = zeros(nGrid, nGrid);
yGrid = zeros(nGrid, nGrid);

for ii=1:nGrid
    for jj=1:nGrid
        xGrid(jj, ii) = grid(ii);
        yGrid(jj, ii) = grid(jj);
        x = grid(ii);
        y = grid(jj);
        z(jj, ii) = (2*pi*sigma(1)*sigma(2)*sqrt(1-rho^2))^(-1)*...
            exp(-(1/(2*(1-rho^2)))*...
            ((x-mu(1))^2/sigma(1)^2 + (y-mu(2))^2/sigma(2)^2-...
            2*rho*(x-mu(1))*(y-mu(2))/(sigma(1)*sigma(2))));
        
    end 
end


time1 = toc;

%%

mesh(xGrid, yGrid, z)


%%

tic;

xGrid = repmat(grid, numel(grid), 1);
yGrid = repmat(grid', 1, numel(grid));

z = (2*pi*sigma(1)*sigma(2)*sqrt(1-rho^2))^(-1)*...
    exp(-(1/(2*(1-rho^2)))*...
    ((xGrid-mu(1)).^2./sigma(1).^2 + (yGrid-mu(2)).^2/sigma(2)^2-...
    2*rho*(xGrid-mu(1)).*(yGrid-mu(2))/(sigma(1)*sigma(2))));

time2 = toc;

%%
mesh(xGrid, yGrid, z)

%%
[xGrid2, yGrid2] = meshgrid(grid, grid);

%%
covVar = [sigma(1)^2 rho*sigma(1)*sigma(2);
    rho*sigma(1)*sigma(2) sigma(2)^2];

z = mvnpdf([xGrid(:) yGrid(:)], mu', covVar);

zMatr = reshape(z, numel(grid), numel(grid));

%%

mesh(xGrid, yGrid, zMatr)

%% strings

s = 'This is string!';

%%

s = ['abc'; 'cba']
%%

s = ['abc'; 'abcd']

%%

c = {[4,5,6,7,9]; 'dkkdkd'; 'sdlkfj'}


%%

stock1.name = 'BMW';
stock1.dates = ['01.01.2001'; '03.02.2001'];
stock1.prices = [112; 114; 111];

%%
stock1.data = []



%%

load patients

patientTable = table(Gender, Age, Smoker, Diastolic);

% subset
subTable = patientTable(:, 1:2);

%
valuesThemself = patientTable{:, 2:4};







