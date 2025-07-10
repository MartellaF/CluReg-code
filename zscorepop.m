function z = zscorepop(x)
if isequal(x,[]), z = []; return; end

% Figure out which dimension sum will work along.
sz = size(x);
dim = find(sz ~= 1, 1);
if isempty(dim), dim = 1; end

% Need to tile the output of mean and std to standardize X.
tile = ones(1,ndims(x)); tile(dim) = sz(dim);

% Compute X's mean and sd, and standardize it.
warn = warning('off','MATLAB:divideByZero');
xbar = repmat(mean(x), tile);
sd = repmat(std(x,1), tile);
warning(warn)
sd(sd==0) = 1; % don't try to scale constant columns
z = (x - xbar) ./ sd;
