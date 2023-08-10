% function [inputsMDL] = load_InputMDL()
% %LOAD_INPUTMDL Load the model description list corresponding to the model
% %inputs. 
% 
% inputsMDL = ...
%         {
%         'thrust'       , 'kN',      'T_{tot}'           ,'Total thrust';...
%         'de_htp'       , 'rad',     '\delta{e}_{htp}'   ,'HTP Deflection';...
%         'da_stick'     , '-',       '\delta{a}_{stick}' ,'Roll Stick Deflection';...
%         'de_stick'     , '-',       '\delta{e}_{stick}' ,'Pitch Stick Deflection';...
%         'dr_stick'     , '-',       '\delta{r}_{stick}' ,'Rudder Stick Deflection';...
%         'w_Wg'         , 'm/s',     '{w}_{Wg}'          ,'Vertical wind speed';...
%         };
% 
% 
% end

function [inputMDL] = load_InputMDL(modelName)%, statesMDL)
%LOAD_inputMDL Load the model description list corresponding to the model inputs.
%   Expected implementation in model: There should be a subsystem named 
%   'Inputs' at the model's top level. 'Inputs' should contain a set of
%   'InputBlocks' (ideally from the library) , each of which is connected to an inport.
%   NB: order of external inports MUST correspond to order of internal
%   inports.

% modelName = 'sim_flexible_unsteady_dev';

inputBlockNames = []; %List of names of types of input blocks. Names must match name of corresponding block in Simulink model/library.
inputBlockIdx = 1;

% Stick inputs
inputBlockNames = [inputBlockNames; {'StickCmd'}];
inputBlockMDL{inputBlockIdx} = ...
        {
        'da_stick'     , '-',       '\delta{a}_{stick}' ,'Roll Stick Deflection';...
        'de_stick'     , '-',       '\delta{e}_{stick}' ,'Pitch Stick Deflection';...
        'dr_stick'     , '-',       '\delta{r}_{stick}' ,'Rudder Stick Deflection';...                
        };

% Basic Flight Controls
inputBlockIdx = inputBlockIdx +1; 
inputBlockNames = [inputBlockNames; {'FltCtrlBasic'}];
inputBlockMDL{inputBlockIdx} = ...
        {
        'throttle'       , '-',      '{\delta}F'           ,'Total thrust';...
        'de_htp'       , 'rad',     '\delta{e}_{htp}'   ,'HTP Deflection';...               
        };

% External wind % TODO: extend to x- and y- directions. 
inputBlockIdx = inputBlockIdx +1; 
inputBlockNames = [inputBlockNames; {'ExtWind'}];
inputBlockMDL{inputBlockIdx} = ...
        {
        'w_Wg'         , 'm/s',     '{w}_{Wg}'          ,'Vertical wind speed';...             
        };
    
% Vertical external wind divided into zones 
inputBlockIdx = inputBlockIdx +1; 
inputBlockNames = [inputBlockNames; {'ExtWindZones'}];
inputBlockMDL{inputBlockIdx} = ...
        {
        'w_Wg_main_wing'   , 'm/s',     '{w}_{Wg,wing}'         ,'Vertical wind speed at foremost main wing section';...    
        'w_Wg_htp'         , 'm/s',     '{w}_{Wg,htp}'          ,'Vertical wind speed at foremost HTP section';...  
        'w_Wg_vtp'         , 'm/s',     '{w}_{Wg,vtp}'          ,'Vertical wind speed at foremost VTP section';...  
        'w_Wg_fuse'        , 'm/s',     '{w}_{Wg,fus}'          ,'Vertical wind speed at foremost fuselage section';...  
        };

% Primary control surface inputs (ailerons not ganged)
inputBlockIdx = inputBlockIdx +1; 
inputBlockNames = [inputBlockNames; {'PrimaryCmd'}];
inputBlockMDL{inputBlockIdx} = ...
        {
        'da_l',	'-',	'\delta_{\xi,l}',   'Left aileron command'; ...
        'da_r',	'-',    '\delta_{\xi,r}',   'Right aileron command'; ...
        'de',	'-',	'\delta_{\eta}',    'Elevator command'; ...
        'dr',	'-',	'\delta_{\zeta}',   'Rudder command'...
        };

% Primary control surface commanded inputs with standard coupling for
% control
inputBlockIdx = inputBlockIdx +1; 
inputBlockNames = [inputBlockNames; {'PrimCtrlCmd'}];
inputBlockMDL{inputBlockIdx} = ...
        {
        'da_opp_cmd',   '-',	'\delta a_{c,opp}',   'Opposed aileron command'; ...
        'da_sym_cmd',	'-',    '\delta a_{c,sym}',   'Symmetric aileron command'; ...
        'de_cmd',       '-',	'\delta e_c',         'Elevator command'; ...
        'dr_cmd',       '-',	'\delta r_c',         'Rudder command'...
        };
    
% Load alleviation devices commands
inputBlockIdx = inputBlockIdx +1; 
inputBlockNames = [inputBlockNames; {'LAD cmd'}];
inputBlockMDL{inputBlockIdx} = ...
        {
        'dLADs',	'-',	'\delta_{LADs}',   'Local load alleviation devices commands' ...
        };
    
%% Construct final MDL:

feval(modelName,[],[],[],'compile');

try
    %Obtain number of outports and line handles
    inportNames = find_system([modelName, '/Inputs'], 'SearchDepth', 1,  'BlockType', 'Inport');
    numInports = size(inportNames,1);
    inportLineHandles = get_param(inportNames, 'LineHandles');
catch
    numInports = 0;
end

inputMDL = []; % Initialize final MDL

for in = 1:numInports
    modelBlockHandle_tmp = get_param(inportLineHandles{in}.Outport, 'dstBlockHandle'); % Find block connected to current inport
    modelBlockName_tmp = get_param(modelBlockHandle_tmp, 'Name'); % Find name of said block
    
    % Compare name of block to all defined MDL blocks 
    inputBlockIdx = find(strcmp(modelBlockName_tmp, inputBlockNames));
    
    if ~isempty(inputBlockIdx)
        if size(inputBlockIdx,1)>1
            error(['Input Model Description List block name "',modelBlockName_tmp{end},'" is not unique. Check output model description list definition.']);
        else
            ph = get_param(modelBlockHandle_tmp,'PortHandles');
            dim = get_param(ph.Inport,'CompiledPortDimensions');
            len = prod(dim);
            if len > size(inputBlockMDL{inputBlockIdx},1)
                for i = 1:len
                    inputBlockMdlAppend = inputBlockMDL{inputBlockIdx};
                    inputBlockMdlAppend{1} = [inputBlockMdlAppend{1},'__',num2str(i)];
                    inputMDL = [inputMDL; inputBlockMdlAppend];
                end
            else
                inputMDL = [inputMDL; inputBlockMDL{inputBlockIdx}]; % Add the matched MDL block to the rest of the final MDL
%             inputBlockNames{inputBlockIdx}=[]; % Remove the inputBlock name from the list so it can't be used twice. 
            end
        end
    else
        error(['"', modelBlockName_tmp,'" does not correspond to any defined output Model Description List block.']);
    end

end

feval(modelName,[],[],[],'term');

end