%% Assistance_fitting
clear all; close all; clc;
addpath(genpath('../'));

%% Subject Data
subject_data = readtable('./logs/subject_data.csv');
pID = "P818";
patient_index = find(strcmp(subject_data.ID, pID));

if isempty(patient_index)
    error("Patient ID was not found!");
    return;
else
    sex = subject_data.sex{patient_index}; % char, 'm' or 'f'
    weight = subject_data.weight(patient_index); % double
    height = subject_data.height(patient_index); % double
end

fname = "./logs/"+pID+"/"+pID+"_IDF.txt";

% Note: Since there is a tab delimiter at the end of each line, an extra
% column of data was being created. These options ignore that column.
% In the future, we should write out the data with comma delimiters (not
% whitespace, especially tabs...) and have newline ends, not delimiter
% ends.
opts = detectImportOptions(fname);
opts.ExtraColumnsRule = 'ignore'; % ignore extra columns created by extra delimiters at the end of lines

data_raw = readtable(fname, opts);
data_raw.Properties.VariableNames = {'Time' 'Elevation' 'Reference' 'Force1' 'Force2' 'Assistance' 'Motor'};

% Hard coded parameters related to data acquisition frequency
clock_frequency = 1000; % Hz
sampling_frequency = 200; % Hz
ts = clock_frequency/sampling_frequency; % ms

%% Force DC Offset
% Force2 tends to have a DC offset, and sometimes it is negative (which
% does not make much sense). For now we just remove this DC offset (though
% in the future we need to take care of modelling cable pre-tension).

% For now I will implement the offset computation naively by just averaging
% some fixed portion of the data. This assumes that the system is at rest
% longer than the window of data that is used to compute the DC offset.
dc_window_len_sec = 5;
dc_window_n = floor(dc_window_len_sec * 1000 / ts);
data_raw.Force1 = data_raw.Force1 - mean(data_raw.Force1(1:dc_window_n));
data_raw.Force2 = data_raw.Force2 - mean(data_raw.Force2(1:dc_window_n));

% Add Force2 into Force1 to see what happens ;)
%data_raw.Force1 = data_raw.Force1 + data_raw.Force2;

%% Fill Missing Data w/ Linear Interpolation
% Generate a new time vector without missing elements.
% This method assumes that we are trying to fill ALL gaps in the data.
time_filled = [data_raw.Time(1):ts:data_raw.Time(end)]';

elevation_filled_ts = resample(timeseries(data_raw.Elevation, data_raw.Time), time_filled, 'linear');
reference_filled_ts = resample(timeseries(data_raw.Reference, data_raw.Time), time_filled, 'linear');
force1_filled_ts = resample(timeseries(data_raw.Force1, data_raw.Time), time_filled, 'linear');
force2_filled_ts = resample(timeseries(data_raw.Force2, data_raw.Time), time_filled, 'linear');
assistance_filled_ts = resample(timeseries(data_raw.Assistance, data_raw.Time), time_filled, 'linear');
motor_filled_ts = resample(timeseries(data_raw.Motor, data_raw.Time), time_filled, 'zoh');

time_filled = time_filled / 1000; % convert to seconds

data_filled = array2table([time_filled elevation_filled_ts.data, reference_filled_ts.data, force1_filled_ts.data, force2_filled_ts.data, assistance_filled_ts.data, motor_filled_ts.data]);
data_filled.Properties.VariableNames = data_raw.Properties.VariableNames; % copy column headers from original table

%% Smoothened Derivative of Elevation
[b,g] = sgolay(5,101); % 101 points with dt = 0.005s corresponds to roughly a 0.5s window
dt = ts/1000;
n_deriv = 1;
dx = zeros(length(data_filled.Elevation),n_deriv+1);
for p = 0:n_deriv
  dx(:,p+1) = conv(data_filled.Elevation, factorial(p)/(-dt)^p * g(:,p+1), 'same');
end

% A plot that shows the elevation angle and derivative estimate
figure; hold on; grid on;
yyaxis left;
plot(data_filled.Time, dx(:,1) ,'b-');
ylabel('Elevation Angle (deg)');

yyaxis right;
plot(data_filled.Time, dx(:,2) ,'m-');
ylim([-20 20])
ylabel('Elevation Angular Speed (deg/s)');

xlabel('Time (s)');

d_th_elev_thresh = 2; % deg/s
ind_slow = find(abs(dx(:,2)) < d_th_elev_thresh);

% This table only contains the data when the arm is moving slowly.
data_slow = data_filled(ind_slow,:);

%% FIT Data
% Visualize with filled data
figure; hold on; grid on;
yyaxis left;
plot(data_filled.Time, data_filled.Force1, 'color', [.8 .6 .8], 'linewidth', 1.1)

plot(data_filled.Time, data_filled.Force2, 'color', [.7 .7 .9], 'linewidth', 1.1)
ylabel('Force (N)');

yyaxis right;
plot(data_filled.Time, data_filled.Reference, 'k')
plot(data_filled.Time, data_filled.Elevation, 'color', [.7 .3 .4], 'linewidth', 1.3, 'linestyle', '-')
ylabel('Elevation Angle (°)');

legend('Force1', 'Force2', 'Reference', 'Elevation', 'location', 'best');
ylim([0 100])
xlabel('Time (s)')
title('Collected Data for Parameter Estimation')

%% Force Correlations
%{
% This scatter plot also does a good job of showing that the two cable
% tension forces are highly correlated (they tend to make a fairly straight
% line when plotted against each other).
figure; hold on; grid on;
plot(data_filled.Force1, data_filled.Force2,'b.');
axis equal;
xlabel('Force1 (N)');
ylabel('Force2 (N)');
title('Scatter Plot of Forces');

[xcorr_forces, lags_forces] = xcorr(data_filled.Force1, data_filled.Force2,'normalized');
lags_forces = lags_forces * ts / 1000; % from milliseconds to seconds
n_downsample = 500; % i.e. every 500 points will result in a spacing of 2.5s with ts = 5 

% Wide xcorr function shows us that the data is highly correlated.
figure; hold on; grid on;
stem(downsample(lags_forces,n_downsample), downsample(xcorr_forces,n_downsample));
xlabel('Signal Lag (s)');
ylabel('Normalized Signal Cross-Correlation');
title('Cable Tension Cross-Correlation Stem Plot');
%}

%% Non-linear curve fitting
% Here we can pick between the raw data, the filled data, and the slow
% data.
%data = data_raw;
%data = data_filled;
data = data_slow;

% Copy data from table
angle = deg2rad(data.Elevation);
force = data.Force1;

% ASSISTIVE FORCE MODEL
% x : [a, b1, b2, c, d1, d2]
F = @(x,theta) (x(1)*sin(theta))  ./  sqrt(1-(x(2)+x(3)*cos(theta+x(4))).^2 ./ (x(5)+x(6)*cos(theta+x(4))));  

% ANTHROPOMETRY
% weights of body segments
if isequal(sex, 'f') % female
    % section lengths
    % l_u, l_f, l_h
    lengths = [0.159 0.152 0.098]*height;
    % segment radii
    % r_u, r_f, r_h
    radii = [0.148 0.094 0.154].*lengths;
    % centre of mass positions
    % c_u, c_f, c_h
    centers = [0.575 0.456 0.343].*lengths;
    % section mass
    % m_u, m_f, m_h
    mass = [0.0255 0.0138 0.0056]*weight;
elseif isequal(sex, 'm') % 'm', male
    % section lengths
    % l_u, l_f, l_h
    lengths = [0.162 0.155 0.108]*height;
    % segment radii
    % r_u, r_f, r_h
    radii = [0.158 0.121 0.184].*lengths;
    % centre of mass positions
    % c_u, c_f, c_h
    centers = [0.577 0.457 0.362].*lengths;
    % section mass
    % m_u, m_f, m_h
    mass = [0.0271 0.0162 0.0061]*weight;
else
    error("Check participant sex!");
    return;
end

% MYOSHIRT GEOMETRY
% horizontal distance between GH joint and Bowden for medial tendon
l_ab0 = 0.02; % m, Jay-value: 0.1
% vertical distance GH joint to Bowden anchor (Anthropometry and extracted from Solidworks)
l_bc0 = 0.023*height + 0.025; % m
% distance between GH joint and bowden anchor
l_ac0 = sqrt(l_ab0^2+l_bc0^2);
% angle between vertical through GH joint and line connecting GH and Bowden anchor
phi0 = atan(l_ab0/l_bc0);
% distance between GH joint and medial anchor point
l_cd0 = vecnorm([0 -lengths(1)+0.050 -radii(1)]); % m, Jay-value: 0.15;
% gravity force corresponding to shoulder torque at 90° elevation
Fg0 = sum((lengths*[0 1 1;0 0 1;0 0 0]+centers).*mass)'*9.81*sin(pi/2) / l_cd0;

% initial parameters
x0 = [Fg0 ...
      l_cd0 ...
      l_ac0 ...
      phi0 ...
      l_ab0^2+l_bc0^2+l_cd0^2 ...
      2*l_cd0*l_ac0];

% % ALTERNATIVE MODEL - DO NOT USE
% % x : [m_arm, k, h, phi]
F_alt = @(x,theta) (x(1).*sin(theta).*9.81) ./ ( x(3).*sqrt( 1-( (x(3)+x(2).*cos(x(4)+theta)).^2 ./ (x(2).^2+x(3).^2+2.*x(2).*x(3).*cos(x(4)+theta)) ) ) );
k0_alt = l_ac0;
h0_alt = l_cd0;
phi0_alt = phi0;
x0_alt = [ l_cd0*Fg0/9.81 ... % sum(mass)
       k0_alt ...
       h0_alt ...
       phi0_alt];

% upper and lower bound for initial values - within 10%
lb = x0*0.9; %used to be 0.9
ub = x0*1.1; %used to be 1.1

lb_alt = x0_alt * 0.9;
ub_alt = x0_alt * 1.1;

% nonlinear lsq fit
%myfit = lsqcurvefit(F,x0,angle,force,lb,ub);
myfit = lsqcurvefit(F,x0,angle,force); % without bounding

%myfit_alt = lsqcurvefit(F_alt,x0_alt,angle,force,lb_alt,ub_alt);
myfit_alt = lsqcurvefit(F_alt,x0_alt,angle,force); % without bounding

fprintf("\n           x0 = {")
fprintf("%.4f, ",x0(1:end-1))
fprintf("%.4f}\n",x0(end))

fprintf("\nfloat param[] = {")
fprintf("%.4f, ",myfit(1:end-1))
fprintf("%.4f};\n",myfit(end))

fprintf("\n           x0 = {")
fprintf("%.4f, ",x0_alt(1:end-1))
fprintf("%.4f}\n",x0_alt(end))

fprintf("\nfloat param[] = {")
fprintf("%.4f, ",myfit_alt(1:end-1))
fprintf("%.4f};\n",myfit_alt(end))
%%
% ANTHROPOMETRY GRAVITY TORQUE MODEL (Fully extended arm model)
T_g = @(theta)sum((lengths*[0 1 1;0 0 1;0 0 0]+centers).*mass)'*9.81*sin(theta);

%th = 0:1:160;
th = linspace(0,180,200);
T_mod = T_g(deg2rad(th));
F_mod = F(x0,deg2rad(th));
F_mod_alt = F_alt(x0_alt,deg2rad(th));

F_fit = F(myfit,deg2rad(th));
F_fit_alt = F_alt(myfit_alt,deg2rad(th));

% FORCE1
figure; hold on; grid on;

% Raw data of angle vs cable tension
scatter(rad2deg(angle),force,3,'markerfacecolor',[.8 .6 .8],'markeredgecolor','none','markerfacealpha',0.5);

% 6-Parameter model used in Mobility study
plot(th,F_mod,'--','color',[.2 .5 .6],'linewidth',1.3);
plot(th,F_fit,'color',[.6 .3 .4],'linewidth',1.3);

% Alternative 4-parameter model (not used in study, retains vertical
% discontinuity)
%plot(th,F_mod_alt,'--','color',[.2 .5 .2],'linewidth',1.3);
%plot(th,F_fit_alt,'color',[.4 .3 .1],'linewidth',1.3);

legend('Data', 'Initial Model', 'Model Fit', 'location', 'best');
%legend('Data', 'Initial Model', 'Initial Alt. Model', 'Model Fit', 'Alt. Model Fit', 'location', 'best');

xlim([0 180]);
ylim([0 200]);
ylabel('Force1 (N)')
xlabel('Elevation (°)')
title('Data Fit')
