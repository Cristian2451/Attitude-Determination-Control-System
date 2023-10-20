% Usmaan 22-5-23
function [T_a, p_a, rho, g1] = atmosphere(y)
% Using NASA atmosphere model for pressure: https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
% And "definition of standard temperatures at different altitudes" for temperature: http://www.ambrsoft.com/CalcPhysics/Altitude/altitude.htm

p0 = 101325; % pa
R = 287;

% gravitational model 
Re = 6378.1e3; % m
g = 9.80665;  % m/s^2
g1 = g .* (Re ./ (Re+y)).^2; % m/s^2

T0 = 288.15; % Temperature at y = 0
L0 = -0.0065; % Rate of change of temperature with height (in K/m)

T1 = 216.65;
p1 = p0.*(T1/T0)^(-g/(R*L0));

T2 = T1;
L2 = 0.001;
p2 = p1*exp(g*(11000 - 20000)/(R*T2));

T3 = 228.65;
L3 = 0.0028;
p3 = p2*(T3/T2)^(-g/(R*L2));

T4 = 270.65;
p4 = p3*(T4/T3)^(-g/(R*L3));

T5 = T4;
L5 = -0.0028;
p5 = p4*exp(g*(47000 - 51000)/(R*T4));

T6 = 214.65;
L6 = -0.002;
p6 = p5*(T6/T5)^(-g/(R*L5));

T7 = 186.95;
p7 = p6*(T7/T6)^(-g/(R*L6));

T8 = T7;
L8 = 0.004;
p8 = p7*exp(g*(84852 - 90000)/(R*T7));

%Work out atmospheric temperature in K
if y < 11000 %Troposphere (0 to 1)
    T_a = T0 + L0*y;
    p_a = p0*(T_a/T0)^(-g/(R*L0));
elseif y >= 11000 && y < 20000 %Tropopause (1 to 2)
    T_a = T1;
    p_a = p1*exp(g*(11000 - y)/(R*T1));
elseif y >= 20000 && y < 32000 %Lower stratosphere (2 to 3)
    T_a = T2 + L2*(y - 20000);
    p_a = p2*(T_a/T2)^(-g/(R*L2));
elseif y >= 32000 && y < 47000 %Upper stratosphere (3 to 4)
    T_a = T3 + L3*(y - 32000);
    p_a = p3*(T_a/T3)^(-g/(R*L3));
elseif y >= 47000 && y < 51000 %Stratopause (4 to 5)
    T_a = T4;
    p_a = p4*exp(g*(47000 - y)/(R*T4));
elseif y >= 51000 && y < 71000 %Lower mesosphere (5 to 6)
    T_a = T5 + L5*(y - 51000);
    p_a = p5*(T_a/T5)^(-g/(R*L5));
elseif y >= 71000 && y < 84852 %Upper mesosphere (6 to 7)
    T_a = T6 + L6*(y - 71000);
    p_a = p6*(T_a/T6)^(-g/(R*L6));
elseif y >= 84852 && y < 90000 %Mesopause (7 to 8)
    T_a = T7;
    p_a = p7*exp(g*(84852 - y)/(R*T7));
else %Thermosphere
    T_a = T8 + L8*(y - 90000);
    p_a = p8*(T_a/T8)^(-g/(R*L8));
end

rho = p_a/(0.2869*T_a*1000);
end