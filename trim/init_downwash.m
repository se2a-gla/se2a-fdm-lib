function [trimPoint] = init_downwash(trimPoint, downwash, gamma_trim)

V = norm( tpGetState(trimPoint,{'u_Kb','v_Kb','w_Kb'}), 2 );

downwash = downwashUnstSchedule( downwash, V );
[A,B,C,D] = tf2ssOcf( downwash.b, downwash.a );%ss(tf(downwash.b, downwash.a));

x_trim = -A \ B*gamma_trim;

trimPoint = tpAddVal( trimPoint, 'state', 'x_downwash__', x_trim );

trimPoint = tpAddVal( trimPoint, 'state', 'delay_downwash', gamma_trim );

end

