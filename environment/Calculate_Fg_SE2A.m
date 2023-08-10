function [F_g] = Calculate_Fg_SE2A(alt)
% Compute flight profile alleviation factor for SE2A.
% Inputs:
% - alt: Altitude (m)
%TODO: change data based on aircraft model?

%Flight profile alleviation factor
Zmo = 11200;    % (m)
MTOW = 64158;  % (kg)
MZFW = 36121+19650;   % (kg)
MLW = 0.9*MTOW;   % (kg) %Rounded up from  A320: MTOW ~ 75t, MLW ~ 65t

F_g = CS25_341_Fg(alt, Zmo, MTOW, MLW, MZFW);

end

