function [dist,indexOnGoalPoints]= map2Points(inputPoints,goalPoints,criterion,numOfPts)
% [dist,indexOnGoalPoints]= map2Points(inputPoints,goalPoints,criterion,numOfPts)
%
% Map from one point cloud to another point cloud, based on the Euclidean
% distance.
% 
% Note this function will give "out of memory" error if the input and goal
% point clouds are too big (>50K).
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

inputPoints = int16(inputPoints);
goalPoints = int16(goalPoints);

N_inputPoints = size(inputPoints,1);
N_goalPoints = size(goalPoints,1);

N_dim_input = size(inputPoints,2);
N_dim_goal = size(goalPoints,2);
% now supports any number of dimensions

if N_dim_input ~= N_dim_goal
    error('You cannot map points that are in spaces of different dimensions!')
end

temp1 = repmat(reshape(inputPoints',1,N_dim_input,N_inputPoints),[N_goalPoints 1 1]);

temp2 = repmat(goalPoints,[1 1 N_inputPoints]);

dist = squeeze(sqrt(sum((temp1 - temp2).^2,2)));

[dist,ind_sortedDist] = sort(dist);

switch criterion
    case 'closest'
        dist = dist(1,:);
        indexOnGoalPoints = ind_sortedDist(1,:);
    case 'farthest'
        dist = dist(end,:);
        indexOnGoalPoints = ind_sortedDist(end,:);
    case 'closer'
        if isempty(numOfPts), error('You want to map to the first XX closest points, please specify XX in the 4th argument.'); end
        if numOfPts>N_goalPoints, error('Number of points exceed size of goal point cloud.'); end
        dist = dist(1:numOfPts,:);
        indexOnGoalPoints = ind_sortedDist(1:numOfPts,:);
    case 'farther'
        if isempty(numOfPts), error('You want to map to the first XX farthest points, please specify XX in the 4th argument.'); end
        if numOfPts>N_goalPoints, error('Number of points exceed size of goal point cloud.'); end
        dist = dist(end-numOfPts+1:end,:);
        indexOnGoalPoints = ind_sortedDist(end-numOfPts+1:end,:);
    otherwise
        error('Please specify either closest, farthest, closer or farther as the mapping criterion.');
end