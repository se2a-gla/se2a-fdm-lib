function [R_b_i, R_b, Q_b] = enginesGetLoads( engines, thrust, eta, pos_ref_c )
% enginesGetLoads compute thrust vectors per engine, total thrust vector
% and total moment w.r.t. instantaneous cg depending on the deformation.

R_b = zeros(3,1);
Q_b = zeros(3,1);

engines_deflection = engines.T_es * eta;

engines_deflection_6xn = reshape( engines_deflection, 6, [] );

engines_displacement = engines_deflection_6xn(1:3,:);
engines_rotation = engines_deflection_6xn(4:6,:);

engines_direction = zeros( 3, engines.num );
for i = 1:engines.num
    engines_direction(:,i) = euler2Dcm( engines_rotation(:,i) )' * engines.direction(:,i);
end

R_b_i = engines_direction * thrust/engines.num;

for i = 1:engines.num
    r_ref = engines.pos_c(:,i) + engines_displacement(:,i) - pos_ref_c;

    [ R_b_i(:,i), Q_b_i ] = forceMomentTransform( R_b_i(:,i), zeros(3,1), r_ref, eye(3) );
    
    R_b = R_b + R_b_i(:,i);
    Q_b = Q_b + Q_b_i;
end

end