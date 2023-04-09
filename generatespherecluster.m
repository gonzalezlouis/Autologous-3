clearvars 

rng shuffle % sets the rng to change with each new instance of calling it 

% Goal: to build an algorithm that plots spheres inside of a cube of
% dimensions L x W x H in which they do not touch each other or the walls,
% output is a list of triplets corresponding to the centers of the spheres
% in each coordinate

n = 30; % number of spheres
r = 10; % cell radius size in um
dim = [1000 1000 1000]; % dimensions of the box in um

lambda = 3; % cell-cell spacing

% Step 1: first we need to create an empty vector of size n x 3
% as a placeholder

EmptyCellVec = zeros(n,3); 

% Step 2: we define a center cell at [1500 1000 50] that we will use as the
% signal measurer in our experiment, this will be the cell with which we
% record the anisotropy from in COMSOL

CenterCell = [500 500 500];

% Step 3: TreatedVec is a 1x3 vector, but in the end we want a 600x3 list
% of triplets. So we need to run this in a while-loop. We also need to
% include the no-overlap condition. This involves checking whether the
% euclidean distance between two spheres is larger than the diameter
% squared of any one of the spheres. If not, then move on to the next one.
% If yes, then append it to our EmptyCellVec (which we will eventually call
% FinalCellVec).

k = 0;
looptime = 1;
diametersquared = (2 * r) ^ 2;
W = dim - 2 * r;

while k < n && looptime < 1e7
    % k is the variable we initialize this loop with and that we index
    TreatedVec = rand(1, 3) .* W + r;
    EuclDist = sum((EmptyCellVec(1:k,:) - TreatedVec).^2, 2);
    a = abs(TreatedVec - CenterCell);
    
    if all(lambda > diametersquared) && all(a >= 5)
        k = k + 1;
        EmptyCellVec(k,:) = TreatedVec;
    end
    EmptyCellVec(1,:) = CenterCell;
    looptime = looptime + 1;
    
    if k > n % display an error message if not enough spheres were found
        error('Too few values found! Only %d were able to be calculated!', ...
            length(EmptyCellVec))
    end
end

% The loop goes as follows: 1) create an index called k that starts at 0
% and a looptime that starts at 1, 2) prepare a "pre-treatment" to apply to
% the random vectors that we generate to enforce non-overlapping and no
% wall clipping, 3) establish a while-loop that starts at 0 and ends at n
% number of cells with a total looptime up to 1e6, 4) calculate the
% Euclidean distance between each new entry of EmptyCellVec and our newly
% generated TreatedVec, 5) calculate the absolute value of the difference
% between the TreatedVec and the CenterCell, 6) run an if-loop that checks
% that every single entry in EuclDist is greater than (diameter)^2 and if
% the difference from the center sphere is >=10 away (this is to make sure
% there is no overlapping with the center sphere), 7) if successful, then
% the TreatedVec will be appended to EmptyCellVec until k = n while also
% appending the CenterCell as one of the entries in the list

% For Line 39, "sum" was used. From Matlab's documentation: S = sum(A,dim)
% returns the sum along dimension dim. For example, if A is a matrix, then
% sum(A,2) is a column vector containing the sum of each row. The Euclidean
% distance goes as d = sqrt((p1-q1)^2 + ....). To simplify and speed up
% coding, we rewrite this expression so we have a d^2 = sum((p1-q1)^2...).
% We need a column vector in order to logically compare each entry with
% diametersquared, which is a scalar. This d^2 is diametersquared in the
% code.


% Step 4: finally we shift the values back so that the CenterCell is
% exactly at [0 0 0] and we rename "EmptyCellVec" to "FinalCellVec" for
% posterity. We then plot this data as surface plots with each dot
% corresponding to a sphere to visualize the results and we save the data
% in .csv format to import into COMSOL for simulating 

FinalCellVec = EmptyCellVec;
xfix = FinalCellVec(:,1) - 500;
yfix = FinalCellVec(:,2) - 500;
zfix = FinalCellVec(:,3) - 500;

FinalCellVec = [xfix,yfix,zfix];

theta = acos(zfix ./ sqrt(xfix.^2 + yfix.^2 + zfix.^2));

clf;
axes('NextPlot', 'add', ...
    'XLim', [-500, 500], 'YLim', [-500, 500], 'ZLim', [-500, 500]);
view(3);
[X, Y, Z] = sphere();
for k = 1:n
    surf(X * r + FinalCellVec(k, 1), Y * r + FinalCellVec(k, 2), ... 
        Z * r + FinalCellVec(k, 3));
end

% writematrix(FinalCellVec,'../dat/cellcoordinatessphere.csv')
