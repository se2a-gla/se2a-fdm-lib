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
            'eta_f' ,  '-',    ''    , 'Filtered structure deflection', [model_name,'/',controller_block_name,'/INDI GLA Controller/filtered measurements/PT2 discrete with saturation/Discrete-Time',newline,'Integrator y'];...
            };
end

end