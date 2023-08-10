function aircraft = aircraftSteady2Unsteady( aircraft )
% only works if the aircraft was initially created with unsteady
% aerodynamics and then converted with aircraftUnsteady2Steady.

field_names = fieldnames(aircraft);

for i = 1:length(field_names)
    current_field_name = field_names{i};
    if contains(current_field_name,'wing')
        aircraft.(current_field_name).config.is_unsteady = 1;
        aircraft.(current_field_name).config.is_circulation_iteration = 0;
    elseif contains(current_field_name,'fuselage')
        aircraft.(current_field_name).config.is_unsteady = 1;
    end
end

end