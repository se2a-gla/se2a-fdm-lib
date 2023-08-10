function simTaskInitBaseWorkspace(aircraft,controller,trimPoint)
% simTaskInitBaseWorkspace assign all required Simulink model parameters to
% base workspace

assignin('base','aircraft',aircraft);
assignin('base','controller',controller);
assignin('base','trimPoint',trimPoint);

ic = tpGenerateIC(trimPoint);
assignin('base','ic',ic);

envir = envirLoadParams('envir_params_default');
pos_ref = loadParams('reference_position_params_se2a');
pos_ref.alt = 0;
assignin('base','envir',envir);
assignin('base','pos_ref',pos_ref);

if ~isfield(trimPoint, 'mdl_names')
    fp.Altitude = -ic.s_Kg(3);
    fp.V_TAS    = norm(ic.V_Kb,2);
else
    fp.Altitude = -tpGetOutput(trimPoint,'z_Kg');
    fp.V_TAS    = tpGetOutput(trimPoint,'V_TAS');
end

% Atmosphere parameters at current flight condition
Atmosphere = Init_Atmosphere_SE2A(fp.Altitude);
Atmosphere.WindConfig.Fg = Calculate_Fg_SE2A(fp.Altitude);
Atmosphere.WindConfig.WindModelType = 'OneMinusCosineGusts';
Atmosphere.WindConfig.U_ds          = -1; % -1 indicates U_ds should be calculated
Atmosphere.WindConfig.x_start       = 1 * fp.V_TAS;

assignin('base','Atmosphere',Atmosphere);

end