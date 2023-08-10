function [outputMDL] = load_OutputMDL(modelName, statesMDL)
%LOAD_OUTPUTMDL Load the model description list corresponding to the model outputs.
%   Expected implementation in model: There should be a subsystem named 
%   'Outputs' at the model's top level. 'Outputs' should contain a set of
%   'OutputBlocks' (ideally from the library , each of which is connected to an outport. 
%   NB: order of external outports MUST correspond to order of internal
%   outports.

% modelName = 'sim_flexible_unsteady_dev';

outputBlockNames = []; %List of names of types of output blocks. Names must match name of corresponding block in Simulink model/library.
outputBlockIdx = 1;

%6DOF Rigid body
outputBlockNames = [outputBlockNames; {'Rigid6DOF_Out'}];
outputBlockMDL{outputBlockIdx} = ...
        {
        % Inertial position in Earth frame:                        
        'x_Kg'           , 'm',     'x_{k,g}'                   ,'X-Position in Earth Frame ';...  
        'y_Kg'           , 'm',     'y_{k,g}'                   ,'Y-Position in Earth Frame';...   
        'z_Kg'           , 'm',     'z_{k,g}'                   ,'Z-Position in Earth Frame';...
        % Euler angles:
        'Phi'            , 'deg',     '\Phi'                    ,'Roll Angle';...                           
        'Theta'          , 'deg',     '\Theta'                  ,'Pitch Angle';...                         
        'Psi'            , 'deg',     '\Psi'                    ,'Heading Angle';... 
        % Inertial velocity in body-fixed frame:
        'u_Kb'           , 'm/s',     'u_{k,b}'                 ,'Inertial Velocity along X-axis in Body-Fixed Frame';...  
        'v_Kb'           , 'm/s',     'v_{k,b}'                 ,'Inertial Velocity along Y-axis in Body-Fixed Frame';...   
        'w_Kb'           , 'm/s',     'w_{k,b}'                 ,'Inertial Velocity along Z-axis in Body-Fixed Frame';...   
        % Inertial velocity in Earth frame:
        'u_Kg'           , 'm/s',     'u_{k,g}'                 ,'Inertial Velocity along X-axis in Earth Frame';...  
        'v_Kg'           , 'm/s',     'v_{k,g}'                 ,'Inertial Velocity along Y-axis in Earth Frame';...   
        'w_Kg'           , 'm/s',     'w_{k,g}'                 ,'Inertial Velocity along Z-axis in Earth Frame';...  
        % Angular rates:
        'p_Kb'           , 'deg/s',   'p_{k,b}'                 ,'Roll Rate';...                        
        'q_Kb'           , 'deg/s',   'q_{k,b}'                 ,'Pitch Rate';...                         
        'r_Kb'           , 'deg/s',   'r_{k,b}'                 ,'Yaw Rate';...                          
        % Inertial accelerations in body-fixed frame:                     
        'u_Kb_dt'        , 'm/s^2',   '\dot{u}_{k,b}'           ,'Inertial Acceleration along X-axis in Body-Fixed Frame';...  
        'v_Kb_dt'        , 'm/s^2',   '\dot{v}_{k,b}'           ,'Inertial Acceleration along Y-axis in Body-Fixed Frame';...   
        'w_Kb_dt'        , 'm/s^2',   '\dot{w}_{k,b}'           ,'Inertial Acceleration along Z-axis in Body-Fixed Frame';...  
        % Angular acccelerations:
        'p_Kb_dt'        , 'deg/s^2', '\dot{p}_{k,b}'           ,'Roll Acceleration';...                        
        'q_Kb_dt'        , 'deg/s^2', '\dot{q}_{k,b}'           ,'Pitch Acceleration';...                         
        'r_Kb_dt'        , 'deg/s^2', '\dot{r}_{k,b}'           ,'Yaw Acceleration';...                   
        };

%Basic air data    
outputBlockIdx = outputBlockIdx +1;  
outputBlockNames = [outputBlockNames; {'AirDataBasic_Out'}];
outputBlockMDL{outputBlockIdx} = ...
        {
        % Inertial position in Earth frame:                        
        'alpha'          , 'deg',       '\alpha'                  ,'Global Angle of Attack';...  
        'alpha_dt'       , 'deg/s',     '\dot{\alpha}'            ,'Global Angle of Attack Rate';...   
        'beta'           , 'deg',       '\beta'                   ,'Global Angle of Sideslip';...
        % Air speeds
        'V_TAS'          , 'm/s',       'V_{TAS}'                 ,'True Air Speed';...
        'V_EAS'          , 'm/s',       'V_{EAS}'                 ,'Equivalent Air Speed';...
        'Mach'           , '-',         'M'                       ,'Mach';...
        % Atmospheric data
        'a_atm'          , 'm/s',       'a_{atm}'                 ,'Speed of Sound';...  
        'p_atm'          , 'Pa',        'p_{atm}'                 ,'Static Atmospheric Pressure';...   
        'T_atm'          , 'K',         'T_{atm}'                 ,'Outside Air Temperature';...
        'rho_atm'        , 'kg/m^3',    '\rho_{atm}'              ,'Air Density';...  
        % Wind speeds
        'u_Wg'           , 'm/s',       'u_{W,g}'                 ,'Wind Speed along X-axis in Earth frame';...   
        'v_Wg'           , 'm/s',       'v_{W,g}'                 ,'Wind Speed along Y-axis in Earth frame';...
        'w_Wg'           , 'm/s',       'w_{W,g}'                 ,'Wind Speed along Z-axis in Earth frame';...
        
        };

%Basic loads   %TODO: synchronize outputs with load_stations struct
outputBlockIdx = outputBlockIdx +1;  
outputBlockNames = [outputBlockNames; {'LoadsBasic_Out'}];
outputBlockMDL{outputBlockIdx} = ...
        {
        % Right wing bending moment                       
        'WR_Root_MX'     , 'Nm',       'M_{x, WR Root}'           ,'Vertical Bending Moment at Right Wing Root ';...  
        'WR_10p_MX'      , 'Nm',       'M_{x, WR 10\%}'           ,'Vertical Bending Moment at 10% Span Right Wing';...
        'WR_20p_MX'      , 'Nm',       'M_{x, WR 20\%}'           ,'Vertical Bending Moment at 20% Span Right Wing';... 
        'WR_30p_MX'      , 'Nm',       'M_{x, WR 30\%}'           ,'Vertical Bending Moment at 30% Span Right Wing';... 
        'WR_40p_MX'      , 'Nm',       'M_{x, WR 40\%}'           ,'Vertical Bending Moment at 40% Span Right Wing';... 
        'WR_50p_MX'      , 'Nm',       'M_{x, WR 50\%}'           ,'Vertical Bending Moment at 50% Span Right Wing';... 
        'WR_60p_MX'      , 'Nm',       'M_{x, WR 60\%}'           ,'Vertical Bending Moment at 60% Span Right Wing';... 
        'WR_70p_MX'      , 'Nm',       'M_{x, WR 70\%}'           ,'Vertical Bending Moment at 70% Span Right Wing';... 
        'WR_80p_MX'      , 'Nm',       'M_{x, WR 80\%}'           ,'Vertical Bending Moment at 80% Span Right Wing';... 
        'WR_90p_MX'      , 'Nm',       'M_{x, WR 90\%}'           ,'Vertical Bending Moment at 90% Span Right Wing';... 
        
        %Left wing bending moment
        'WL_Root_MX'     , 'Nm',       'M_{x, WL Root}'           ,'Vertical Bending Moment at Left Wing Root ';...  
        'WL_10p_MX'      , 'Nm',       'M_{x, WL 10\%}'           ,'Vertical Bending Moment at 10% Span Left Wing';...
        'WL_20p_MX'      , 'Nm',       'M_{x, WL 20\%}'           ,'Vertical Bending Moment at 20% Span Left Wing';... 
        'WL_30p_MX'      , 'Nm',       'M_{x, WL 30\%}'           ,'Vertical Bending Moment at 30% Span Left Wing';... 
        'WL_40p_MX'      , 'Nm',       'M_{x, WL 40\%}'           ,'Vertical Bending Moment at 40% Span Left Wing';... 
        'WL_50p_MX'      , 'Nm',       'M_{x, WL 50\%}'           ,'Vertical Bending Moment at 50% Span Left Wing';... 
        'WL_60p_MX'      , 'Nm',       'M_{x, WL 60\%}'           ,'Vertical Bending Moment at 60% Span Left Wing';... 
        'WL_70p_MX'      , 'Nm',       'M_{x, WL 70\%}'           ,'Vertical Bending Moment at 70% Span Left Wing';... 
        'WL_80p_MX'      , 'Nm',       'M_{x, WL 80\%}'           ,'Vertical Bending Moment at 80% Span Left Wing';... 
        'WL_90p_MX'      , 'Nm',       'M_{x, WL 90\%}'           ,'Vertical Bending Moment at 90% Span Left Wing';...  
        };

%Flex. structure generalized coordinates
outputBlockIdx = outputBlockIdx +1;  
outputBlockNames = [outputBlockNames; {'FlexStruct_Out'}];

eta_idx_tmp = find(startsWith(statesMDL.Continuous(:,1), 'eta__'));
eta_dt_idx_tmp = find(startsWith(statesMDL.Continuous(:,1), 'eta_dt__'));
eta_dt2_tmp = statesMDL.Continuous(eta_idx_tmp,1:4);
eta_dt2_tmp(:,1) = generate_numbered_str_array('eta_dt2__', 1:length(eta_idx_tmp));
eta_dt2_tmp(:,2) = strcat(eta_dt2_tmp(:,2), '/s^2');
eta_dt2_tmp(:,3) = generate_numbered_str_array('\ddot{\eta}_', 1:length(eta_idx_tmp));
eta_dt2_tmp(:,4) = repmat({'Flex. mode acceleration'}, size(eta_dt2_tmp(:,4)));

outputBlockMDL{outputBlockIdx} = ...
        [statesMDL.Continuous(eta_idx_tmp,1:4);...
         statesMDL.Continuous(eta_dt_idx_tmp,1:4);...
         eta_dt2_tmp...
         ];
clear eta_idx_tmp eta_dt_idx_tmp eta_dt2_tmp


%% Construct final MDL:

try
    %Obtain number of outports and line handles
    outportNames = find_system([modelName, '/Outputs'], 'SearchDepth', 1,  'BlockType', 'Outport');
    numOutports = size(outportNames,1);

    outportLineHandles = get_param(outportNames, 'LineHandles');
catch
    numOutports = 0;
end

outputMDL = []; % Initialize final MDL

for in = 1:numOutports
    modelBlockHandle_tmp = get_param(outportLineHandles{in}.Inport, 'srcBlockHandle'); % Find block connected to current outport
    modelBlockName_tmp = get_param(modelBlockHandle_tmp, 'Name'); % Find name of said block
    
    % Compare name of block to all defined MDL blocks 
    outputBlockIdx = find(strcmp(modelBlockName_tmp, outputBlockNames));
    
    if ~isempty(outputBlockIdx)
        if size(outputBlockIdx,1)>1
            error(['Output Model Description List block name "',modelBlockName_tmp{end},'" is not unique. Check output model description list definition.']);
        else
            outputMDL = [outputMDL; outputBlockMDL{outputBlockIdx}]; % Add the matched MDL block to the rest of the final MDL
%             outputBlockNames{outputBlockIdx}=[]; % Remove the outputBlock name from the list so it can't be used twice. 

        end
    else
        error(['"', modelBlockName_tmp,'" does not correspond to any defined output Model Description List block.']);
    end

end



end