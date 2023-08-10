function trim_reqts = tpCreateReqtsRigid( altitude, V_EAS )


listOfTasks = [];

% Col 1: Type, Col 2: Variable, Col 3: Values
% Col 1 choices: Y_Reqts, Xdot_Reqts, XCont_Var, XDisc_Var, U_Var
sweepData  = {'Y_Reqts',     'z_Kg',      -altitude  
              'Y_Reqts',     'V_EAS',     V_EAS %'Y_Reqts',     'u_Kb',      0
              'Y_Reqts',     'w_Kg',      0
              'Xdot_Reqts',  'u_Kb',      0
              'Xdot_Reqts',  'w_Kb',      0
              'Xdot_Reqts',  'q_Kb',      0
              'Xdot_Reqts',  'dHTP',      0
              'Xdot_Reqts',  'dHTP_dt',      0
              'Xdot_Reqts',  'thrust',      0
              'XCont_Var',  'z__g',      [0]
              'XCont_Var',  'u_Kb',      150
              'XCont_Var',  'w_Kb',      0
              'XCont_Var',  'qt_2',      0
              'XCont_Var',  'dHTP',      0
              'XCont_Var',  'dHTP_dt',      0
              'XCont_Var',  'thrust',      1000
              'U_Var',      'throttle',      0.5
              'U_Var',      'de_htp',      0};

sweepType = 1; % 1 = cross, 2 = parallel 

% Assemble trim cases

trim_reqts = [listOfTasks, Assemble_TrimSweep(sweepData, sweepType)];
% 

end