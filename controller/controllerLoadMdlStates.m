function controller_states = controllerLoadMdlStates(model_name, controller_block_name )


switch controller_block_name
    case 0
        controller_states.Continuous = [];
        controller_states.Discrete = [];
    case 'INDI Controller'
        controller_states.Continuous = [];
        controller_states.Discrete = ...
            { ...
            'eta_d' ,  '-',    ''    , 'Desired structure deflection', [model_name,'/',controller_block_name,'/INDI GLA Controller/setpoints/PT2 discrete with saturation/Discrete-Time',newline,'Integrator y'];...
            'eta_f' ,  '-',    ''    , 'Filtered structure deflection', [model_name,'/',controller_block_name,'/INDI GLA Controller/filtered measurements/Structural deflection filter/Discrete-Time',newline,'Integrator y'];...
            'a_z_Kb_f' ,  '-',    ''    , 'Filtered vertical acceleration', [model_name,'/',controller_block_name,'/INDI GLA Controller/filtered measurements/Vertical acceleration filter/Discrete-Time',newline,'Integrator y'];...
            'v_f', '-', '', 'Filtered velocity', [model_name, '/', controller_block_name, '/Airplane Pitch Rate Flight Mode with INDI/Velocity filter/Discrete-Time',newline,'Integrator y'];...
            'q_f', '-', '', 'Filtered pitch rate', [model_name, '/', controller_block_name, '/Airplane Pitch Rate Flight Mode with INDI/Pitch rate filter/Discrete-Time',newline,'Integrator y'];...
            };
end

end