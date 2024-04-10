function [G,G_local] = indiCeLadVar( cel, rho, V_A, Ma )
% Scheduling of LAD control effectiveness

% Disclamer:
%   SPDX-License-Identifier: GPL-3.0-only
% 
%   Copyright (C) 2023 Yannic Beyer
%   Copyright (C) 2023 TU Braunschweig, Institute of Flight Guidance
% *************************************************************************

V_A = max(V_A,1);


beta_Ma_flap = sqrtReal(1-Ma.^2);
cla = 2*pi./beta_Ma_flap;
F = airfoilFlapEffectiveness(0.15);
delta_dt=0;
delta = 1;
delta_qs = F.F_10.*delta/pi + F.F_11.*delta_dt.*cel.ci./(4*pi*V_A);
dcLdu = cla*delta_qs * 0.5236;


% [1], eq. (51)
delta_qs_M = -(F.F_4+F.F_10)./(2*pi.*beta_Ma_flap).*delta ...
    -((2*F.F_1-2.*F.F_8-(2*F.e+1).*F.F_4+F.F_11)./(8*pi*beta_Ma_flap)) .* cel.ci./V_A.* delta_dt;
dcmdu = pi/beta_Ma_flap*delta_qs_M;

% dcmdu(:) = 0;

coeff2force = rho/2*(V_A.*cel.cos_sweep).^2.*cel.si;
u2force = diag(dcLdu.*coeff2force);
coeff2moment = coeff2force.*cel.ci;
u2moment = diag(dcmdu.*coeff2moment);

G_local = u2force * cel.cntrlim;

cel_acc_z = sum(cel.spanwise_mapping*G_local,1)/cel.m;

cel_lift = cel.daccdl ...
    * u2force ...
    * cel.cntrlim;
% cel_lift(:) = mean(cel_lift);

cel_moment = cel.daccdm ...
    * u2moment ...
    * cel.cntrlim;

% control effectiveness matrix
% G = cel_lift + cel_moment;
G = [cel_acc_z; cel_lift + cel_moment*0];
% G(3,:) = -200*ones(size(G(3,:)));
% G3 = cel.T_sc(6+21,4:4:end).*diag(u2moment)';
% G(3,:) = G3([1:19,22:end]);
% G(3,:) = cel.T_sc(6+21,4:4:end)*u2moment*cel.spanwise_mapping*cel.cntrlim;

end