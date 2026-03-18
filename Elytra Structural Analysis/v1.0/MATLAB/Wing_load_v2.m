clc;
clear;

% all units are SI (m, N, N*m, Pa, etc.)
b = 2.9; % span length
lambda = 1; % taper ratio
C_r = 0.331; % root chord
L = 168.56; % total lift
x_cl_y = 0.136; % approximate wing centerline height from bottom of fuselage
y_fus = 0.152; % approximate fuselage diameter (for purpose of estimation, fuselage is not circular)
t = 0.0025; % fuselage skin thickness
sigma_T_ult_PLA = 11.49*10^6; % fuselage skin material ultimate tensile stress (for 3D printerd PLA using weakest direction)
E_PLA = 1.171*10^9; % PLA elastic modulus in weakest direction
FoS = 1.5; % factor of safety (14 CFR part 23.2230 since part 107 doesn't have a FoS requirement)
load_factor = 2.5; % maximum operational load factor


% Calculate wing lift distribution, moment at wing root, fuselage skin
% thickness, and graph the lift distribution
% -------------------------------------------------------------------------
% The following two functions aren't directly used but are combined in
% function L_bar_yf
% trapezoidal lift as function of span
L_tf = @(y) (2*L)/(b*(1+lambda))*(1-2*(y/b)*(1-lambda));

% elliptical lift as function of span
L_ef = @(y) (4*L)/(pi*b)*sqrt(1-(2*y/b).^2);

function [L_bar_yf, M_bar_y2] = wingLiftMoment(L_cruise, b, lambda, n)
    % average lift distribution (Schrenk's method)
    L = L_cruise*n;
    L_bar_yf = @(y) 0.5.*((2*L)/(b*(1+lambda))*(1-2*(y/b)*(1-lambda))+(4*L)/(pi*b)*sqrt(1-(2*y/b).^2));
    M_bar_y2f = @(y) 0.5.*((2*L)/(b*(1+lambda))*(1-2*(y/b)*(1-lambda))+(4*L)/(pi*b)*sqrt(1-(2*y/b).^2)).*y;
    M_bar_y2 = integral(M_bar_y2f, 0, b/2);
end

% Normal Load
[L_bar_yf_norm, M_bar_y2_norm] = wingLiftMoment(L,b,lambda,1);
[L_bar_yf_op, M_bar_y2_op] = wingLiftMoment(L,b,lambda,load_factor);
[L_bar_yf_ult, M_bar_y2_ult] = wingLiftMoment(L,b,lambda,load_factor*FoS);

L_bar_y = integral(L_bar_yf_norm, -b/2, b/2);
if isapprox(L_bar_y, L)
    disp('Lift distribution passes lift-sum check')
else
    disp('Lift distribution does not pass lift-sum check')
end

% Reaction moment at wing root
% disp(['Reaction Moment at Wing Root: ', num2str(M_bar_y2), ' N*m'])

function [sigma_T_factored, FoSBool] = skinTens(R, M_bar_y2, t, sigma_T_bound, n)
    % Fuselage skin tensile stress
    sigma_T = 2*M_bar_y2/(pi*R^2*t);
    sigma_T_factored=sigma_T;
    % check if factor of safety is met
    FoSBool = (sigma_T_factored < sigma_T_bound);
end

R = y_fus/2;

% Normal load
[sigma_T_norm, FoSBool_skin_norm] = skinTens(R, M_bar_y2_norm, t, sigma_T_ult_PLA, 1);

% Maximum operational load
[sigma_T_op, FoSBool_skin_op] = skinTens(R, M_bar_y2_op, t, sigma_T_ult_PLA, load_factor);

% Ultimate load (max op * FoS)
[sigma_T_ult, FoSBool_skin_ult] = skinTens(R, M_bar_y2_ult, t, sigma_T_ult_PLA, load_factor*FoS);

disp('---------------------------------')
disp('Fuselage skin tensile stress:')
disp(['Normal operation load: ',num2str(sigma_T_norm*10^-6),' MPa'])
disp(['Meets requirements: ', string(FoSBool_skin_norm)])
disp(['Max operation load: ',num2str(sigma_T_op*10^-6),' MPa'])
disp(['Meets requirements: ', string(FoSBool_skin_op)])
disp(['Ultimate load FoS: ',num2str(sigma_T_ult*10^-6),' MPa'])
disp(['Meets requirements: ', string(FoSBool_skin_ult)])
disp('---------------------------------')

% Calculate factor of safety with current outer-wall thickness settings
MoS_skin = sigma_T_ult_PLA/sigma_T_ult-1;
disp(['Calculated worst-case skin thickness Margin of Safety: ', num2str(MoS_skin)])

% Minimum fuselage skin thickness given 1.5 factor of safety
t_min = 2*M_bar_y2_ult/(pi*sigma_T_ult_PLA*R^2);
t_min_mm = t_min*1000;
disp(['Minimum fuselage skin thickness at ultimate load and FoS ',num2str(FoS),': ',num2str(t_min_mm),' mm'])

% toggle figure because it was annoying me
figyes = 0;
if figyes == 1
    figure
    fplot(L_bar_yf_norm, [-b/2, b/2])
    title('Lift distribution')
    xlabel('Span Position (m)')
    ylabel('Lift (N)')
    ylim([0; 80])
end
disp('---------------------------------')
% -------------------------------------------------------------------------

sparDouble = 0;
if sparDouble == 0
    % single 2020 extrusion
    h = 0.02; % spar height, extruded AL 2020 is square
    I_AL = 6.826*10^-9; % AL 2020 moment of inertia (same in both directions)
    A_spar = 0.0001591;
elseif sparDouble == 1
    % two 2020 extrusions stacked and bonded tightly
    h = 0.02*2;
    A_spar = 0.0001591; % total area of the new spar (from 80/20 website)
    I_AL = 2*(6.826*10^-9+A_spar*(0.01)^2);
else
    % 2040 extrusion
    h = 0.04;
    I_AL = 4.5357*10^-8;
end
sigma_T_yield_AL = 172.37*10^6; % from 80/20 website
sigma_T_ult_AL = 194*10^6; % based on the reduced yield strength of the 2020 compared to normal 6063-T6 Al
E_AL = 69*10^9; % from 80/20 website

% Calculate spar tensile and torsional stress
% We don't know the bending moment proportion absorbed by the spar, so
% assume it absorbs all of it
% -------------------------------------------------------------------------
%sigma_T_AL_max = M_bar_y2_norm*(h/2)/I_AL;
%sigma_T_AL_max_MPA = sigma_T_AL_max/10^6;
%disp(['Maximum Tensile Stress in Spar: ', num2str(sigma_T_AL_max_MPA), ' MPa']);

% calculate FoS of spar tensile stress
%n_des_spar_tensile = sigma_T_yield_AL/sigma_T_AL_max;
%disp(['Calculated AL 2020 spar FoS: ',num2str(n_des_spar_tensile)])

function [sigma_T_AL_max, FoSBool, tau_2020] = sparStress(L, M_bar_y2,h,I_AL,sigma_T_bound,n,A_spar)
    sigma_T_AL_max = M_bar_y2*(h/2)/I_AL;
    FoSBool = (sigma_T_AL_max < sigma_T_bound);
    tau_2020 = (L*n/2)*3/(2*A_spar);
end

[sigma_T_AL_max_norm, FoSBool_spar_norm, tau_2020_norm] = sparStress(L, M_bar_y2_norm, h, I_AL, sigma_T_yield_AL, 1, A_spar);
[sigma_T_AL_max_op, FoSBool_spar_op, tau_2020_op] = sparStress(L, M_bar_y2_op, h, I_AL, sigma_T_yield_AL, load_factor, A_spar);
[sigma_T_AL_max_ult, FoSBool_spar_ult, tau_2020_ult] = sparStress(L, M_bar_y2_ult, h, I_AL, sigma_T_ult_AL, load_factor*FoS, A_spar);

disp('Wing spar tensile stress:')
disp(['Normal operation load: ',num2str(sigma_T_AL_max_norm*10^-6),' MPa'])
disp(['Meets requirements: ', string(FoSBool_spar_norm)])
disp(['Max operation load: ',num2str(sigma_T_AL_max_op*10^-6),' MPa'])
disp(['Meets requirements: ', string(FoSBool_spar_op)])
disp(['Ultimate load FoS: ',num2str(sigma_T_AL_max_ult*10^-6),' MPa'])
disp(['Meets requirements: ', string(FoSBool_spar_ult)])
disp('---------------------------------')

disp(['Spar yield margin of safety: ',num2str(sigma_T_yield_AL/sigma_T_AL_max_op-1)])
disp(['Spar ultimate margin of safety: ',num2str(sigma_T_ult_AL/sigma_T_AL_max_ult-1)])
disp('---------------------------------')
% calculate section modulus
% S_2020 = I_AL / (h / 2);
% S_allow = M_bar_y2_op/(sigma_T_yield_AL/load_factor);
% disp(['Is section modulus of selected spar acceptable? ',num2str(S_2020),'>',num2str(S_allow)])
% if S_2020>S_allow
%     disp('----> yes')
% else
%     disp('----> no')
% end

% calculate shear stress in spar
% tau_2020 = (L_bar_y/2)*3/(2*h^2);
disp(['Maximum shear stress of the spar: '])
disp(['Normal operation load: ',num2str(tau_2020_norm*10^-6),' MPa'])
disp(['Max operation load: ',num2str(tau_2020_op*10^-6),' MPa'])
disp(['Ultimate load FoS: ',num2str(tau_2020_ult*10^-6),' MPa'])
disp('---------------------------------')

% calculate deflection at tip
function [delta_spar_tip] = sparDeflect(L_bar_yf, E_AL, I_AL, b)
    delta_spar_tip_f = @(y) L_bar_yf(y)/(E_AL*I_AL).*(b/2-y).^2/2;
    delta_spar_tip = integral(delta_spar_tip_f,0,b/2);
end

[delta_spar_tip_norm] = sparDeflect(L_bar_yf_norm, E_AL, I_AL, b);
[delta_spar_tip_op] = sparDeflect(L_bar_yf_op, E_AL, I_AL, b);
[delta_spar_tip_ult] = sparDeflect(L_bar_yf_ult, E_AL, I_AL, b);

disp('Deflection at tip of spar: ')
disp(['Normal operation load: ', num2str(delta_spar_tip_norm), ' m']);
disp(['----> As percentage of semi-span: ', num2str(delta_spar_tip_norm/(b/2)*100), '%', newline]);
disp(['Max operation load: ', num2str(delta_spar_tip_op), ' m']);
disp(['----> As percentage of semi-span: ', num2str(delta_spar_tip_op/(b/2)*100), '%', newline]);
disp(['Ultimate load FoS: ', num2str(delta_spar_tip_ult), ' m']);
disp(['----> As percentage of semi-span: ', num2str(delta_spar_tip_ult/(b/2)*100), '%']);
disp('---------------------------------')
