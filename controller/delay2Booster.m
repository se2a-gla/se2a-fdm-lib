function [boost_servo,boost_aero] = delay2Booster( T_delay, omega_sens, dtc, atc )

T_delay_sens = 2/omega_sens;
T_delay_act = T_delay - T_delay_sens;
boost_servo = dtc / (T_delay_act/2);
boost_aero = min(atc) / (T_delay_act/2);

end