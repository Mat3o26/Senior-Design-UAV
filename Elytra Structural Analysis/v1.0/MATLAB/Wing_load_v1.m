clc;
clear;

% all units are SI (m, N, N*m, Pa, etc.)
b = 2.9; % span length
lambda = 1; % taper ratio
C_r = 0.331; % root chord
L = 230; % total lift
x_cl_y = 0.136; % approximate wing centerline height from bottom of fuselage
y_fus = 0.152; % approximate fuselage diameter (for purpose of estimation, fuselage is not circular)
t = 0.0025; % fuselage skin thickness
sigma_T_ult_PLA = 11.49*10^6; % fuselage skin material ultimate tensile stress (for 3D printerd PLA using weakest direction)
n_req = 1.5; % factor of safety


% Calculate wing lift distribution, moment at wing root, fuselage skin
% thickness, and graph the lift distribution
% -------------------------------------------------------------------------
% The following two functions aren't directly used but are combined in
% function L_bar_yf

% trapezoidal lift as function of span
L_tf = @(y) (2*L)/(b*(1+lambda))*(1-2*(y/b)*(1-lambda));

% elliptical lift as function of span
L_ef = @(y) (4*L)/(pi*b)*sqrt(1-(2*y/b).^2);

% average lift distribution (Schrenk's method)
L_bar_yf = @(y) 0.5.*((2*L)/(b*(1+lambda))*(1-2*(y/b)*(1-lambda))+(4*L)/(pi*b)*sqrt(1-(2*y/b).^2));
L_bar_y = integral(L_bar_yf, -b/2, b/2);
if isapprox(L_bar_y, L)
    disp('Lift distribution passes lift-sum check')
else
    disp('Lift distribution does not pass lift-sum check')
end

% Reaction moment at wing root
M_bar_y2f = @(y) 0.5.*((2*L)/(b*(1+lambda))*(1-2*(y/b)*(1-lambda))+(4*L)/(pi*b)*sqrt(1-(2*y/b).^2)).*y;
M_bar_y2 = integral(M_bar_y2f, 0, b/2);
disp(['Reaction Moment at Wing Root: ', num2str(M_bar_y2), ' N*m'])

% Fuselage skin tensile stress
R = y_fus/2;
% I = pi/2*R^3*t;
sigma_T = 2*M_bar_y2/(pi*R^2*t);
sigma_T_MPa = sigma_T/10^6;
disp(['Fuselage skin tensile stress: ',num2str(sigma_T_MPa),' MPa'])

% Calculate factor of safety with current outer-wall thickness settings
n_des_skin = sigma_T_ult_PLA/sigma_T;
disp(['Calculated Worst-Case Factor of Safety: ', num2str(n_des_skin)])

% Minimum fuselage skin thickness given 1.5 factor of safety
t_min = 2*M_bar_y2*n_req/(pi*sigma_T_ult_PLA*R^2);
t_min_mm = t_min*1000;
disp(['Minimum fuselage skin thickness with FoS ',num2str(n_req),': ',num2str(t_min_mm),' mm'])

% toggle figure because it was annoying me
figyes = 0;
if figyes == 1
    figure
    fplot(L_bar_yf, [-b/2, b/2])
    title('Lift distribution')
    xlabel('Span Position (m)')
    ylabel('Lift (N)')
    ylim([0; 100])
end
disp('---------------------------------')
% -------------------------------------------------------------------------


h = 0.02; % spar height, extruded AL 2020 is square
I_AL = 6.826*10^-9; % AL 2020 moment of inertia (same in both directions)
sigma_T_yield_AL = 172.37*10^6; % from 80/20 website
E_AL = 69*10^9; % from 80/20 website

% Calculate spar tensile and torsional stress
% We don't know the bending moment proportion absorbed by the spar, so
% assume it absorbs all of it
% -------------------------------------------------------------------------
sigma_T_AL_max = M_bar_y2*(h/2)/I_AL;
sigma_T_AL_max_MPA = sigma_T_AL_max/10^6;
disp(['Maximum Tensile Stress in Spar: ', num2str(sigma_T_AL_max_MPA), ' MPa']);

% calculate section modulus
% 
S_2020 = I_AL / (h / 2);
S_allow = M_bar_y2/(sigma_T_yield_AL/n_req);
disp(['Is section modulus of selected spar acceptable? ',num2str(S_2020),'>',num2str(S_allow)])
if S_2020>S_allow
    disp('----> yes')
else
    disp('----> no')
end

% calculate FoS of spar tensile stress
n_des_spar_tensile = sigma_T_yield_AL/sigma_T_AL_max;
disp(['Calculated AL 2020 spar FoS: ',num2str(n_des_spar_tensile)])

% calculate shear stress in spar
tau_2020 = (L_bar_y/2)*3/(2*h^2);
disp(['maximum shear stress of the spar: ',num2str(tau_2020*10^-6),' MPa'])

% calculate deflection at tip
delta_spar_tip_f = @(y) L_bar_yf(y)/(E_AL*I_AL).*(b/2-y).^2/2;
delta_spar_tip = integral(delta_spar_tip_f,0,b/2);
disp(['Deflection at tip of spar: ', num2str(delta_spar_tip), ' m']);
disp(['Deflection at tip of spar as percentage of semi-span: ', num2str(delta_spar_tip/(b/2)*100), '%']);