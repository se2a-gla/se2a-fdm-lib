function [ic] = tpGenerateIC(trimPoint_in)
%TPGENERATEIC Generates struct 'ic' from a trimPoint

%Default:
if ~isfield(trimPoint_in, 'mdl_names')
    ic.omega_Kb = [ 0; 0; 0 ];
    ic.q_bg     = [ 1; 0; 0; 0 ];
    ic.V_Kb     = [ 150; 0; 0 ];
    ic.s_Kg     = [ 0; 0; 0 ];
else
    ic.omega_Kb = [tpGetState(trimPoint_in, 'p_Kb'); tpGetState(trimPoint_in, 'q_Kb'); tpGetState(trimPoint_in, 'r_Kb')];
    ic.q_bg     = [tpGetState(trimPoint_in, 'qt_0'); tpGetState(trimPoint_in, 'qt_1'); tpGetState(trimPoint_in, 'qt_2'); tpGetState(trimPoint_in, 'qt_3')];
    ic.V_Kb     = [tpGetState(trimPoint_in, 'u_Kb'); tpGetState(trimPoint_in, 'v_Kb'); tpGetState(trimPoint_in, 'w_Kb')];
    ic.s_Kg     = [tpGetState(trimPoint_in, 'x__g'); tpGetState(trimPoint_in, 'y__g'); tpGetState(trimPoint_in, 'z__g')];
end

end

