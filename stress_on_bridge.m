function bridge_1_driver
% This is the driver fxn

% Display adjacency matrix for nos = 0, 1, 2
for nos = 0:3
    [adj, xc, yc, len] = build_basic_bridge(nos);
    disp(adj)
end

% Plot spy plots for nos = 1, 2, 3
for nos = 1:3
    [adj, xc, yc, len] = build_basic_bridge(nos);
    plot_spy(adj, nos)
    figure
end

% Plot UN_DEFORMED (basic) bridges for nos = 1, 2, 3
car_weight = 0;
for nos = 1:3
    [adj, xc, yc, len] = build_basic_bridge(nos);
    plot_bridge(nos,xc,yc,car_weight)
    figure
end

% Plot DEFORMED (complete) bridges for nos = 1, 2, 3 and weight = 0.01
car_weight = 0.01;
for nos = 1:3
    [adj, dx, dy] = build_complete_bridge(nos, car_weight);
    plot_bridge(nos,dx,dy,car_weight)
    figure
end

% Plot DEFORMED (complete) bridges for nos = 1, 2, 3 and weight = 0.05
car_weight = 0.05;
for nos = 1:3
    [adj, dx, dy] = build_complete_bridge(nos, car_weight);
    plot_bridge(nos,dx,dy,car_weight)
    figure
end

% Plot graphs for Questions #1-4
% Question #1 graphs:
car_weight = 0.01;
for nos = 3:5
    [adj, dx, dy] = build_complete_bridge(nos, car_weight);
    plot_bridge(nos,dx,dy,car_weight)
    figure 
end

car_weight = 0.05;
for nos = 1:3
    [adj, dx, dy] = build_complete_bridge(nos, car_weight);
    plot_bridge(nos,dx,dy,car_weight)
    
    % To prevent the creation of a blank final figure
    if nos ~= 3
        figure
    end
end


% Question Answers:

% #1: The maximum "safe" length (in my opinion) for a bridge to support a
% car weighing 0.01 units is 4 sections. This is because a 5 section bridge
% produces a curve in the bridge that is too steep for cars to brake
% safely, and increases the likelihood of a crash significantly. Cars going
% down the slope into the dip of the bridge would likely accelerate too
% quickly and make it difficult to brake, whereas a bridge with 4 sections
% or less has less of a slope and therefore a more safe experience
% traveling the bridge.

% #2: Using the same logic, an 18-wheeler truck would likely have
% difficulty stopping on a bridge with 2 or more sections due to the
% steepness of the bridge's curve, meaning that the maximum "safe" length
% for an 18-wheeler truck should be 1 section. 

% #3: Generally, the more non-zero points there are on the spy plot, the greater the
% deformation (not including the influence of the weight of the
% car). There is a matrix of points that is 5 rows long and 8 columns wide
% that repeats itself the same number as the number of sections in the
% bridge (figures 1-3). 

% #4: As mentioned in #3, according to the spy plot, the repeating matrix of points repeats itself the
% same number of times as the number of sections in the bridge. It shifts
% itself to the right several units and down by 5 units, and starts from
% the 3rd column. Additionally, the bridge deforms significantly 
% when the number of sections of the bridge
% also increases. The length of each fiber also increases more
% significantly. Evidence for this is in figures 16-18 (weight = 0.05), 
% 10-15 (weight = 0.01) The first set of figures and the second set of
% figures show that the both the length of each fiber and the deformation
% of each bridge increases with the number of sections in the bridge.
end

function [adj, dx, dy] = build_complete_bridge(nos,car_weight)
% Builds a bridge under the weight of cars.
% Inputs:
% nos - number of sections (square sections with fibers structured like a
% cross)
% car_weight - the weight of the cars on the bridge
% Outputs:
% adj - the adjacency matrix of a bridge with nos sections
% dx - a matrix with "number of fibers" rows and 2 columns that for each
% row, has the x-coordinate of a first node (connected by the fiber of that
% row) and the x-coordinate of the last node
% dy - a matrix with "number of fibers" rows and 2 columns that for each
% row, has the y-coordinate of a first node (connected by the fiber of that
% row) and the y-coordinate of the last node


% Introduce adj, xc, yc, len, force variables
[adj, xc, yc, len] = build_basic_bridge(nos);
force = zeros(4 + 4 * nos, 1);
force(2:4:end) = -car_weight;

[dx,dy,work,X,Y] = deform_basic_bridge(nos,adj,xc,yc,len,force);
end

function [adj, xc, yc, len] = build_basic_bridge(nos)
% Builds a bridge without any weights.
% Inputs:
% nos - number of sections (square sections with fibers structured like a
% cross)
% Outputs:
% adj - the adjacency matrix of a bridge with nos sections
% xc - a matrix with "number of fibers" rows and 2 columns that for each
% row, has the x-coordinate of a first node (connected by the fiber of that
% row) and the x-coordinate of the last node
% yc - a matrix with "number of fibers" rows and 2 columns that for each
% row, has the y-coordinate of a first node (connected by the fiber of that
% row) and the y-coordinate of the last node
% len - a matrix with "number of fibers" rows and 1 column such that the
% element in each row is the length of the fiber of the row (fiber) it is
% in

% calculate some helpful numbers
num_nodes = 2 + 2 * nos;
num_edges = 5 + 5 * nos;
s = 1/sqrt(2);
max_x = nos + 2;
% initialize the return values
adj = zeros(num_edges, 2*num_nodes);
xc = zeros(num_edges, 2);
yc = zeros(num_edges, 2);
len = ones(num_edges, 1);
% build the left side of bridge
adj(1, 1) = 1;
adj(2, [3 4]) = [s s];
xc(1, :) = [0 1];
xc(2, :) = [0 1];
yc(1, :) = [0 0];
yc(2, :) = [0 1];
len(2) = 1/s;
% build the middle of bridge
for i = 1:nos
    % 8 x 5 matrix
    M = [0 -1 0 1 0 0 0 0;
        0 0 -s s s -s 0 0;
        0 0 -1 0 0 0 1 0;
        -s -s 0 0 0 0 s s;
        -1 0 0 0 1 0 0 0];
    % 5 * i - 2 is the start for each section
    adj(5 * i - 2: 5 * i + 2, 4 * i - 3: 4 * i + 4) = M;
    xc(5 * i - 2,:) = [i i];
    xc(5 * i - 1: 5 * i + 2, :) = [i i + 1; i i + 1; i i + 1; i i + 1];
    yc(5 * i - 2: 5 * i + 2, :) = [0 1 ; 1 0; 1 1; 0 1; 0 0];
    len(5 * i - 1) = 1/s;
    len(5 * i + 1) = 1/s;
end 

% build the right side of bridge< ... >;
adj(end - 2, end - 2) = -1;
adj(end - 2, end) = 1;
adj(end - 1, [(end - 1) end]) = [-s s];
adj(end, end - 3) = -1;
xc(end - 2, :) = [max_x - 1 max_x - 1];
xc(end - 1, :) = [max_x - 1 max_x];
xc(end, :) = [max_x - 1 max_x];
yc(end - 2, :) = [0 1];
yc(end - 1, :) = [1 0];
yc(end, :) = [0 0];
len(end - 1) = 1/s;
end

function [dx,dy,work,X,Y] = deform_basic_bridge(nos,adj,xc,yc,len,force)
% A helper function that determines the displacement of a bridge based on
% the position of the original nodes of a bridge
% Inputs:
% nos - number of sections (square sections with fibers structured like a
% cross)
% adj - the adjacency matrix of a bridge with nos sections
% xc - a matrix with "number of fibers" rows and 2 columns that for each
% row, has the x-coordinate of a first node (connected by the fiber of that
% row) and the x-coordinate of the last node
% yc - a matrix with "number of fibers" rows and 2 columns that for each
% row, has the y-coordinate of a first node (connected by the fiber of that
% row) and the y-coordinate of the last node
% len - a matrix with "number of fibers" rows and 1 column such that the
% element in each row is the length of the fiber of the row (fiber) it is
% in
% force - a matrix with (4 + 4 * nos) rows and 1 column based on the weight
% of the cars that represents the force on the bridge
% Output:
% dx - a matrix with the same dimensions as xc whose elements represent the
% changed coordinates of the bridge nodes after applying a force
% dy - a matrix with the same dimensions as yc whose elements represent the
% changed coordinates of the bridge nodes after applying a force
% work - a scalar that measures the way the bridge system as a whole has
% deformed and how much work the fibers are doing to prevent this
% deformation
% X - matrices that associate each node with its displacement in the
% horizontal direction
% Y - matrices that associate each node in the bridge with its displacement
% in the vertical direction



% calculate some helpful numbers
stiffness = adj'*diag(1./len)*adj;
displacements = stiffness\force;
work = displacements'*force;
X = displacements(1:2:end);
Y = displacements(2:2:end);
% initialize the return values
dx = zeros(size(xc));
dy = zeros(size(yc));
% deform the left side of bridge
dx(1,:) = xc(1,:) + [0 X(1)];
dx(2,:) = xc(2,:) + [0 X(2)];
dy(1,:) = yc(1,:) + [0 Y(1)];
dy(2,:) = yc(2,:) + [0 Y(2)];
% deform the middle of bridge
for i = 1:nos
    dx(5*i-2,:) = xc(5*i-2,:) + [X(2*i-1) X(2*i)];
    dx(5*i-1,:) = xc(5*i-1,:) + [X(2*i) X(2*i+1)];
    dx(5*i,:) = xc(5*i,:) + [X(2*i) X(2*i+2)];
    dx(5*i+1,:) = xc(5*i+1,:) + [X(2*i-1) X(2*i+2)];
    dx(5*i+2,:) = xc(5*i+2,:) + [X(2*i-1) X(2*i+1)];
    
    dy(5*i-2,:) = yc(5*i-2,:) + [Y(2*i-1) Y(2*i)];
    dy(5*i-1,:) = yc(5*i-1,:) + [Y(2*i) Y(2*i+1)];
    dy(5*i,:) = yc(5*i,:) + [Y(2*i) Y(2*i+2)];
    dy(5*i+1,:) = yc(5*i+1,:) + [Y(2*i-1) Y(2*i+2)];
    dy(5*i+2,:) = yc(5*i+2,:) + [Y(2*i-1) Y(2*i+1)];
end
% deform the right side of bridge
dx(end-2,:) = xc(end-2,:) + [X(end-1) X(end)];
dx(end-1,:) = xc(end-1,:) + [X(end) 0];
dx(end,:) = xc(end,:) + [X(end-1) 0];

dy(end-2,:) = yc(end-2,:) + [Y(end-1) Y(end)];
dy(end-1,:) = yc(end-1,:) + [Y(end) 0];
dy(end,:) = yc(end,:) + [Y(end-1) 0];
end

function plot_bridge(nos,xc,yc,car_weight)
% Plots a bridge given xc and yc (matrices of the nodes' coordinates) and
% other inputs used to detail the graph title.
% Inputs: 
% nos - number of sections (square sections with fibers structured like a
% cross)
% xc - a matrix with "number of fibers" rows and 2 columns that for each
% row, has the x-coordinate of a first node (connected by the fiber of that
% row) and the x-coordinate of the last node
% yc - a matrix with "number of fibers" rows and 2 columns that for each
% row, has the y-coordinate of a first node (connected by the fiber of that
% row) and the y-coordinate of the last node
% car_weight - the weight of the cars on the bridge

hold on
for i = 1:length(xc)
    plot(xc(i,:),yc(i,:),"-o")
end
title("Plot of Bridge with " + num2str(nos) + "Sections and " + num2str(car_weight) + " Car Weight")
xlabel("X-Coordinates of Bridge Nodes")
ylabel("Y-Coordinates of Bridge Nodes")
hold off
end

function plot_spy(adj, nos)
% Plots the sparsity pattern of the adjacency matrix of a bridge with "nos" sections. 
% Inputs: 
% adj - the adjacency matrix of a bridge with nos sections
% nos - number of sections (square sections with fibers structured like a
% cross)
hold on
spy(adj, "o")
title("Spy Plot of Adjacency Matrix of Bridge With " + num2str(nos) + " Sections")
xlabel("X- and Y- Coordinates of Bridge Nodes")
ylabel("Bridge Fibers")
hold off
end