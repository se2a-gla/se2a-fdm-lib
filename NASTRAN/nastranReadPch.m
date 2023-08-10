%******************************************************************************************************
%*                                                                                                    *
%*        (c) 2020      Deutsches Zentrum f�r Luft- und Raumfahrt e.V.                                *
%*                               in der Helmholtz-Gemeinschaft                                        *
%*                                                                                                    *
%*                           Institute of Flight Systems, Braunschweig                                *
%*                                                                                                    *
%******************************************************************************************************
%*                                                                                                    *
%*                                   readPchFromNastran                                               *
%*                                                                                                    *
%******************************************************************************************************
%*                                                                                                    *
%* Type     : Function                                                                                *
%*                                                                                                    *
%* Purpose  : Extract NASTRAN .pch files containing stiffness and mass matrices into struct           * 
%*                                                                                                    *
%* Version  : 1.1                                                                                     *
%*                                                                                                    *
%******************************************************************************************************
%*                                                                                                    *
%* Remarks  : -                                                                                       *
%*                                                                                                    *                                                                                              *
%******************************************************************************************************
%* Author               *       Date            *       Description                                   *
%******************************************************************************************************
%* Davide Cavaliere           August 2020              davide.cavaliere@dlr.de                        *
%* Yannic Beyer               September 2020           y.beyer@tu-bs.de                               *
%******************************************************************************************************

function nastran_pch = nastranReadPch( filename )

disp('nastranReadPch: Starting to read mass and stiffness from .pch file...')

fileID = fopen(filename);

cIdxs{1} = [1:8];
cIdxs{2} = [9:16];
cIdxs{3} = [17:24];
cIdxs{4} = [25:32];
cIdxs{5} = [33:40];
cIdxs{6} = [41:56];
cIdxs{7} = [57:72];

fileRow = 1;
tline = fgetl(fileID);

sysMat = cell(2,1);
sysNodes = [];

while tline ~= -1  
    
    %Remove whitespace and create a cell vector containing contents of row
    for in = 1:6
        colCont{in} = strtrim(tline(cIdxs{in}));
    end
    if length(tline) > 56
        colCont{7} = strtrim(tline(cIdxs{7}));
    end
    
    %Save each row in a single cell matrix
    pchCellMat{fileRow} = colCont;
    
    % Identify which type of line it is
    if strcmp(colCont{1}, 'DMIG')
        lineCase = 1; % Beginning of new matrix
    elseif strcmp(colCont{1}, 'DMIG*')
        lineCase = 2; % Beginning of new node
    else
        lineCase = 3; % Contents of a node
    end
    
    switch lineCase
        case 1 % New matrix is declared
            if strcmp(colCont{2}, 'KAAX')
                matType = 1;  %Stiffness matrix
            elseif strcmp(colCont{2}, 'MAAX')
                matType = 2;  %Mass matrix
            end
            
            if strcmp(colCont{4}, '6')
                matSym = 1;  %Symmetric matrix
            end
            
            nCol    = str2num(colCont{7});
            nNodes  = nCol / 6;
            
            if matSym == 1
                sysMat{matType,1} = zeros(nCol,nCol);
            end
            
            currNode = 0;
            currNodeCol = 0;
            
            nodeCount = 1;
            
        case 2 %Next Node element
            lastNode = currNode;
            lastNodeCol = currNodeCol;
            
            currNode = str2num(colCont{5});
            currNodeCol = str2num(colCont{6}); 
            
            
            if lastNode ~= currNode % New node
                if ~any(currNode == sysNodes) % Node already in list of nodes?
                    sysNodes(nodeCount) = currNode;
                end
                nodeCount = nodeCount+1;
            end
            
            activeCol = (find(currNode==sysNodes)-1)*6 + currNodeCol;
            
        case 3 %Values for the current node element
            tgtNode     = str2num(colCont{3});
            tgtNodeRow  = str2num(colCont{5});
            
            activeRow = (find(tgtNode == sysNodes)-1)*6 + tgtNodeRow;
            
            cellValue = str2num(colCont{6});
            
            %Assign matrix values to matrices:
            sysMat{matType,1}(activeRow, activeCol) = cellValue;
            
            if (activeRow ~= activeCol) && (matSym==1) % Matrix is indeed symmetrical and we're not on the diagonal
                sysMat{matType,1}(activeCol, activeRow) = cellValue;
            end
            
    end
    tline = fgetl(fileID);
    fileRow = fileRow + 1;
end

KAAX = sysMat{1,1};       % Stiffness matrix
MAAX = sysMat{2,1};       % Mass matrix

nodesList   = sysNodes;   % All the node IDs present in this file
numNodes    = nNodes;     % number of nodes
numCol      = nCol;       % # of rows/columns in each matrix
allContents = pchCellMat; % Total contents of the pch file, row by row.

nastran_pch.K         = KAAX;
nastran_pch.M         = MAAX;
nastran_pch.nodesList = nodesList;

fclose(fileID);

disp('nastranReadPch: Finished!')

end
