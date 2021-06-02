%% Assistance_fitting
clear all; close all; clc;
addpath(genpath('../'));

%% Subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pID = "P946";
height = 1.70; %m
weight = 66; % kg
gender = 'male'; % 'female','male'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pID = "P490";
%height = 1.63; %m
%weight = 58; % kg
%gender = 'female'; % 'female','male'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = "./logs/"+pID+"/"+pID+"_IDF.txt";

% Note: Since there is a tab delimiter at the end of each line, an extra
% column of data was being created. These options ignore that column.
opts = detectImportOptions(fname);
opts.ExtraColumnsRule = 'ignore';   %ignore extra columns created by extra delimiters at the end of lines

data = readtable(fname, opts);
data.Properties.VariableNames = {'Time' 'Elevation' 'Reference' 'Force1' 'Force2' 'Assistance' 'Motor'};

%% Extract data
% motor stops when sequence finishes >> extract data from end of array
clock_frequency = 1000; % Hz
sampling_frequency = 200; % Hz
ts = clock_frequency/sampling_frequency;
elevationPeriod = 5000/ts; % ms

%% FIT Data
figure; hold on; grid on;
yyaxis left;
plot(data.Time,data.Force1,'color',[.8 .6 .8],'linewidth',1.1)

plot(data.Time,data.Force2,'color',[.7 .7 .9],'linewidth',1.1)
ylabel('Force (N)');

yyaxis right;
plot(data.Time,data.Reference,'k')
plot(data.Time,data.Elevation,'color',[.7 .3 .4],'linewidth',1.3,'linestyle','-')
ylabel('Elevation Angle (°)');

legend('Force1','Force2','Reference','Elevation','location','best');
ylim([0 100])
xlabel('Time (s)')
title('Collected Data for Parameter Estimation')

%% Elevation Speed


%% Force Correlations
%{
% This scatter plot also does a good job of showing that the two cable
% tension forces are highly correlated (they tend to make a fairly straight
% line when plotted against each other).
figure; hold on; grid on;
plot(data.Force1, data.Force2,'b.');
axis equal;
xlabel('Force1 (N)');
ylabel('Force2 (N)');
title('Scatter Plot of Forces');

[xcorr_forces, lags_forces] = xcorr(data.Force1, data.Force2,'normalized');
lags_forces = lags_forces * ts / 1000; % from milliseconds to seconds
n_downsample = 500; % i.e. every 500 points will result in a spacing of 2.5s with ts = 5 

% Wide xcorr function shows us that the data is highly correlated.
figure; hold on; grid on;
stem(downsample(lags_forces,n_downsample), downsample(xcorr_forces,n_downsample));
%}

%% Non-linear curve fitting

% USE ONLY REST PERIOD DATA
% angle = deg2rad([rest30.Elevation;rest60.Elevation;rest90.Elevation]);
% reference = [rest30.Reference;rest60.Reference;rest90.Reference];
% force = [rest30.Force1;rest60.Force1;rest90.Force1];

% USE ALL DATA - DO THIS FOR NOW
angle = deg2rad(data.Elevation);
reference = data.Reference;
force = data.Force1;

% ASSISTIVE FORCE MODEL
% x : [a, b1, b2, c, d1, d2]
F = @(x,theta) (x(1)*sin(theta))  ./  sqrt(1-(x(2)+x(3)*cos(theta+x(4))).^2 ./ (x(5)+x(6)*cos(theta+x(4))));  

% >>> INITIALIZATION <<<

% ANTHROPOMETRY
% weights of body segments
if isequal(gender, 'female') % female
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
else % 'male', male
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
myfit = lsqcurvefit(F,x0,angle,force,lb,ub);
%myfit = lsqcurvefit(F,x0,angle,force); % without bounding

myfit_alt = lsqcurvefit(F_alt,x0_alt,angle,force,lb_alt,ub_alt);
%myfit_alt = lsqcurvefit(F_alt,x0_alt,angle,force); % without bounding

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
% SAVE TO FILE
%{
filesave = 'params/'+pID+'_params.txt';
fileID = fopen(filesave,'a+');
fprintf(fileID,'= %s ===============\n',pID);
fprintf(fileID,'weight: %.1f kg\n',weight);
fprintf(fileID,'height: %.1f m\n',height);
fprintf(fileID,'gender: %s\n',gender);
fprintf(fileID,'= INITIALIZATION =====\n');
fprintf(fileID,"x0 = {");
fprintf(fileID,"%.4f, ",x0(1:end-1));
fprintf(fileID,"%.4f}\n",x0(end));
fprintf(fileID,'= PARAMETER FIT  =====\n');
fprintf(fileID,"float param[] = {");
fprintf(fileID,"%.4f, ",myfit(1:end-1));
fprintf(fileID,"%.4f};\n\n",myfit(end));
fclose(fileID);
%%
%}
% ANTHROPOMETRY GRAVITY TORQUE MODEL
T_g = @(theta)sum((lengths*[0 1 1;0 0 1;0 0 0]+centers).*mass)'*9.81*sin(theta);

cutVal = find(round(data.Reference)==90);
%th = 0:1:160;
th = linspace(0,180,200);
T_mod = T_g(deg2rad(th));
F_mod = F(x0,deg2rad(th));
F_mod_alt = F_alt(x0_alt,deg2rad(th));

F_fit = F(myfit,deg2rad(th));
F_fit_alt = F_alt(myfit_alt,deg2rad(th));

figure; hold on; grid on;
% FORCE
%plot(data.Elevation(first:last),data.Force1(first:last),'.','color',[.8 .6 .8]);
scatter(rad2deg(angle),force,3,'markerfacecolor',[.8 .6 .8],'markeredgecolor','none','markerfacealpha',0.5);
plot(th,F_mod,'--','color',[.2 .5 .6],'linewidth',1.3);
%plot(th,F_mod_alt,'--','color',[.2 .5 .2],'linewidth',1.3);

plot(th,F_fit,'color',[.6 .3 .4],'linewidth',1.3);
%plot(th,F_fit_alt,'color',[.4 .3 .1],'linewidth',1.3);
%text(20,F_mod(60),'Initialization','color',[.2 .5 .6])
%text(30,F_fit(60),'Fitted','color',[.6 .3 .4])
%text(40,force(ceil(length(force)/8)),'Data','color',[.8 .6 .8])

legend('Data', 'Initial Model', 'Model Fit', 'location', 'best');
%legend('Data', 'Initial Model', 'Initial Alt. Model', 'Model Fit', 'Alt. Model Fit', 'location', 'best');

xlim([0 180]);
ylim([0 200]);
ylabel('Force (N))')
xlabel('Elevation (°)')
title('Data fit')

% TORQUE
% plot(data.Elevation(first:last),data.Force1(first:last).*l_cd0.*sind(data.Elevation(first:last)),'.','color',[.8 .6 .8]);
% plot(th,F_mod.*l_cd0.*sind(th),'color',[.6 .3 .4]);
% plot(th,F_fit.*l_cd0.*sind(th),'color',[.2 .5 .6]);
% plot(th,T_mod,'color','k');

