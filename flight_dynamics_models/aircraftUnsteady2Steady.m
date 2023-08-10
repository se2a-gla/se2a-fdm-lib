function aircraft = aircraftUnsteady2Steady( aircraft )

field_names = fieldnames(aircraft);

for i = 1:length(field_names)
    current_field_name = field_names{i};
    if contains(current_field_name,'wing')
        aircraft.(current_field_name).config.is_unsteady = 0;
        aircraft.(current_field_name).config.is_circulation_iteration = 1;
    elseif contains(current_field_name,'fuselage')
        aircraft.(current_field_name).config.is_unsteady = 0;
    end
end

end