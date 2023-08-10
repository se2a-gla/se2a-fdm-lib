function aircraft = aircraftSe2aUpdateMach( aircraft, Ma )

aircraft.wing_main = wingUpdateMach( aircraft.wing_main, Ma );
aircraft.wing_htp = wingUpdateMach( aircraft.wing_htp, Ma );
aircraft.wing_vtp = wingUpdateMach( aircraft.wing_vtp, Ma );

end