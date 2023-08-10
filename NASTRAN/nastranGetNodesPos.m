function nodesPos = nastranGetNodesPos( xPos, yPos, zPos, nodeIDs, nodesList )
% nastranGetNodesPos returns the position of all that are used for the structure
%   dynamics in the .pch file. Therefore, the position of the GRID is
%   required including the IDs as well as the nodes list.
% 
% Inputs:
%	xPos        position vector in x direction from the GRID
%	yPos        position vector in y direction from the GRID
%	zPos        position vector in z direction from the GRID
% 	nodeIDs     ID vector corresponding to the xPos, yPos, zPos vectors
%   nodeList    vector of node IDs that are used in the .pch file
% 
% Outputs:
%   nodesPos    struct with fields:
%               x   node position vector in x direction
%               y   node position vector in y direction
%               z   node position vector in z direction
% 
% See also: nastranReadGrid
% 
% Authors:
%   Yannic Beyer
% 
%   Copyright 2020 TU Braunschweig
% *************************************************************************

    disp('nastranGetNodesPos: Starting to convert grid to node positions...')

    numNodes = length(nodesList);
    
    % initialize variables
    xPosNode = zeros(numNodes,1);
    yPosNode = zeros(numNodes,1);
    zPosNode = zeros(numNodes,1);
    numIsempty = 0;
    numDuplicate = 0;
    
    % search for equivalent IDs
    for i = 1:numNodes
        idx = find(nodeIDs == nodesList(i));
        % if there are duplicate GRID nodes, take the first one and
        % increment counter
        if length(idx) > 1
            idx = idx(1);
            numDuplicate = numDuplicate + 1;
        end
        % if node is found in the GRID nodes, save its position, else:
        % increment counter
        if ~isempty(idx)
            xPosNode(i) = xPos(idx);
            yPosNode(i) = yPos(idx);
            zPosNode(i) = zPos(idx);
        else
            numIsempty = numIsempty + 1;
        end
    end
    
    disp(['nastranGetNodesPos: Number of nodes found: ',num2str(numNodes)]);
    disp(['nastranGetNodesPos: Number of nodes NOT found: ',num2str(numIsempty)]);
    disp(['nastranGetNodesPos: Number of DUPLICATE nodes: ',num2str(numDuplicate)]);
    
    nodesPos.x = xPosNode;
    nodesPos.y = yPosNode;
    nodesPos.z = zPosNode;
    
    disp('nastranGetNodesPos: Finished!')
    
end