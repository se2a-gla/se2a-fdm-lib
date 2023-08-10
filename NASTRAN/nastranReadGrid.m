%******************************************************************************************************
%*                                                                                                    *
%*        (c) 2020      Deutsches Zentrum f�r Luft- und Raumfahrt e.V.                                *
%*                               in der Helmholtz-Gemeinschaft                                        *
%*                                                                                                    *
%*                           Institute of Flight Systems, Braunschweig                                *
%*                                                                                                    *
%******************************************************************************************************
%*                                                                                                    *
%*                                  readGridFromNastran                                               *
%*                                                                                                    *
%******************************************************************************************************
%*                                                                                                    *
%* Type     : Function                                                                                *
%*                                                                                                    *
%* Purpose  : Extract NASTRAN GRID files containing node positions into grid struct                   * 
%*                                                                                                    *
%* Version  : 1.2                                                                                     *
%*                                                                                                    *
%******************************************************************************************************
%*                                                                                                    *
%* Remarks  : This script assumes all GRID files are in a subdirectory (specified by the folder name).*
%*                                                                                                    *                                                                                              *
%******************************************************************************************************
%* Author               *       Date            *       Description                                   *
%******************************************************************************************************
%* Davide Cavaliere           August 2020              davide.cavaliere@dlr.de                        *
%* Yannic Beyer               September 2020           y.beyer@tu-bs.de                               *
%******************************************************************************************************

function nastran_grid = nastranReadGrid( folder_name )

disp('nastranReadGrid: Starting to read grid from .grid files...')

% Find folder name in path
matlab_path = path;
len_name = length(folder_name);
idx_on_path = strfind( matlab_path, folder_name );
dir_len = 1;
for i = 1:length(idx_on_path)
    current_symbol = matlab_path(idx_on_path(i) + len_name);
    if ~ ( strcmp(current_symbol,';') || strcmp(current_symbol,':') )
        if i == length(idx_on_path)
            error('folder not detected')
        end
        continue;
    end
    while true
        current_symbol = matlab_path(idx_on_path(i) - dir_len);
        if ~ ( strcmp(current_symbol,';') || strcmp(current_symbol,':') )
            dir_len = dir_len + 1;
        else
            dir_len = dir_len - 1;
            break;
        end
    end
    break;
end
idx_begin = idx_on_path(i) - dir_len;
idx_end = idx_on_path(i) + len_name - 1;
% This is the path of the specified folder name
folder_path = matlab_path( idx_begin:idx_end );


%Find all GRID
gridFiles = dir(folder_path);

firstColIdxs = [1:8];
nodeIDIdxs = [9:16];
xPosIdxs = [25:32];
yPosIdxs = [33:40];
zPosIdxs = [41:48];

nodeIDs = [];
xPos = [];
yPos = [];
zPos = [];

for in = 1:length(gridFiles)
    
    name = gridFiles(in).name;
    
    % only consider files that
    if strcmp(name(1),'.') || ~contains(name,'GRID')
        continue;
    end
    
    filename = [gridFiles(in).folder, '/', name];

    fileID = fopen(filename);

    fileRow = 1;
    tline = fgetl(fileID);

    while tline ~= -1  % -1 is EOF
        
        if length(tline) >= 48 % All useful lines have at least 48 characters
            firstCol = strtrim(tline(firstColIdxs));
        else
            firstCol = '';
        end
        
        %Remove whitespace and fill in node coordinates
        if strcmp(firstCol, 'GRID') % All lines with nodes start with GRID
            nodeIDs = [nodeIDs; str2double(strtrim(tline(nodeIDIdxs)))];
            xPos = [xPos; str2doubleExp(strtrim(tline(xPosIdxs)))];
            yPos = [yPos; str2doubleExp(strtrim(tline(yPosIdxs)))];
            zPos = [zPos; str2doubleExp(strtrim(tline(zPosIdxs)))];
        end

        tline = fgetl(fileID);
        fileRow = fileRow + 1;
    end
    fclose(fileID);
end

%Sort nodes in ascending order
[nodeIDs, sortIdxs] = sort(nodeIDs);
xPos = xPos(sortIdxs);
yPos = yPos(sortIdxs);
zPos = zPos(sortIdxs);

nastran_grid.xPos = xPos;
nastran_grid.yPos = yPos;
nastran_grid.zPos = zPos;
nastran_grid.nodeIDs = nodeIDs;

    function double_ = str2doubleExp( string_ )
        idx_sign = 0;
        if contains(string_(2:end),'+')
            idx_sign = strfind( string_(2:end), '+' ) + 1;
        elseif contains(string_(2:end),'-')
            idx_sign = strfind( string_(2:end), '-' ) + 1;
        end
        if idx_sign > 0
            string_2 = [string_(1:idx_sign-1),'e',string_(idx_sign:end)];
        else
            string_2 = string_;
        end
        double_ = str2double(string_2);
    end

disp('nastranReadGrid: Finished!')


end
