% CEE 246 Homework 2 - Response History & Spectrum Analysis
%
% Ray Abbiatici
% Version 1.0/RJA/20-Jan-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

clear all
clc

%% Fundamental Period, Frequency and Spectral Period range
T1 = 1.037;       % [sec]
w1 = (2*pi()/T1); % [Hz]
T = xlsread('_SearchResults.csv',1,'A169:A279'); % [sec]
g = 386.4;        % [in/sec2]

%% Call pSA Spectrum Data from PEER Database
% RSN864
RSN864AX1 = g*xlsread('_SearchResults.csv',1,'B169:B279'); % [in/sec2]
RSN864AX2 = g*xlsread('_SearchResults.csv',1,'C169:C279'); % [in/sec2]

% RSN1616
RSN1616AX1 = g*xlsread('_SearchResults.csv',1,'K169:K279'); % [in/sec2]
RSN1616AX2 = g*xlsread('_SearchResults.csv',1,'L169:L279'); % [in/sec2]

% RSN6915
RSN6915AX1 = g*xlsread('_SearchResults.csv',1,'Z169:Z279');   % [in/sec2]
RSN6915AX2 = g*xlsread('_SearchResults.csv',1,'AA169:AA279'); % [in/sec2]

%% Compute Pseudo Displacement Spectra from Pseudo Acceleration Spectra

for i = 1:length(T)
    w(i) = 2*pi()/T(i);
end

for i = 1:length(w)
    % RSN864
    RSN864DX1(i) = (1/w(i)^2).*RSN864AX1(i);    % [in]
    RSN864DX2(i) = (1/w(i)^2).*RSN864AX2(i);    % [in]
    
    % RSN1616
    RSN1616DX1(i) = (1/w(i)^2).*RSN1616AX1(i);  % [in]
    RSN1616DX2(i) = (1/w(i)^2).*RSN1616AX2(i);  % [in]
    
    % RSN6915
    RSN6915DX1(i) = (1/w(i)^2).*RSN6915AX1(i);  % [in]
    RSN6915DX2(i) = (1/w(i)^2).*RSN6915AX2(i);  % [in]
end

%% Plot Acceleration & Displacement Spectra
figure('Name','Acceleration & Displacement Spectra')
subplot(1,2,1)
plot(T,abs(RSN864AX1))
hold on
plot(T,abs(RSN864AX2))
hold on
plot(T,abs(RSN1616AX1))
hold on
plot(T,abs(RSN1616AX2))
hold on
plot(T,abs(RSN6915AX1))
hold on
plot(T,abs(RSN6915AX2))
xlim([0 3])
xline(T1)
xlabel('T (sec)')
ylabel('Acceleration (in/sec2)')
legend('RSN864X1','RSN864X2','RSN1616X1','RSN1616X2','RSN6915X1', ...
    'RSN6915X2')

subplot(1,2,2)
plot(T,RSN864DX1)
hold on
plot(T,RSN864DX2)
hold on
plot(T,RSN1616DX1)
hold on
plot(T,RSN1616DX2)
hold on
plot(T,RSN6915DX1)
hold on
plot(T,RSN6915DX2)
xlim([0 3])
xline(T1)
xlabel('T (sec)')
ylabel('Displacement (in)')
legend('RSN864X1','RSN864X2','RSN1616X1','RSN1616X2','RSN6915X1', ...
    'RSN6915X2')

disp('The Controlling Ground Motions Are RSN86490, RSN6915W,RSN1616N')

%% Compute the Elastic Response of the System

z = 0.05;        % [% Damping]
M1 = 7.508;      % [k-sec^2/in]
W = (M1/0.67)*g; % [k]
hi = 12.*[0; 15; 27; 39; 51; 63; 75; 87; 99]; % [in]
phi = [0.000; 0.054; 0.164; 0.323; 0.518; 0.739; 0.977; 1.225; 1.476];

%   RSN864X2
run RSN864_LANDERS_JOS090_a;
m = M1;
ag = column(ag);
p = m*ag*g;
c = 2*z*M1*w1;
k = M1*w1^2;
tRSN864 = [0:dt:dt*(length(ag)-1)]';

% For Linear Response
gamma = 1/2;
beta = 1/4;
[aRSN864,vRSN864,uRSN864] = NewmarkIntegrator(gamma,beta,m,c,k,p,dt);
aRSN864 = abs(aRSN864-p/m);

% Base Shear V(t) Normalized by W
VRSN864n = M1*(w1^2)/W*uRSN864;     % [k/k]

% Single DOF and Roof Level Displacements normalized by SDOF and Roof hts
d = zeros(length(hi),length(uRSN864));

for i = 1:length(hi)
    d(i,:) = phi(i).*uRSN864';
end

dtsn = uRSN864 ./ (76.112*12);        % [in/in]
dtrn = d(9,:)' ./ (99*12);            % [in/in]

% Story Drift
drift = zeros([length(hi) length(uRSN864)]);
maxdrift = zeros([length(hi) 1]);

for i = 2:length(hi)
    drift(i,:) = abs(d(i,:) - d (i-1,:))/(hi(i)-hi(i-1))*100; % [%]
    maxdrift(i) = max(drift(i,:));
end 

fprintf('\nBelow are the results of the Response History Analysis:\n')
fprintf('\nPeak Response Quantities for RSN864090 are:\n')
disp('Normalized Base Shear [k/k] =')
disp(max(abs(VRSN864n)))
disp('Base Shear [k] =')
disp(max(abs(VRSN864n)*W))
disp('Overturning Moment [ft-k] =')
disp(max(abs(VRSN864n)*W)*76.112)
disp('Normalized SDOF Displacement [in/in] =')
disp(max(abs(dtsn)))
disp('SDOF Displacement [in] =')
disp(max(abs(dtsn)*76.112*12))
disp('Normalized Roof Displacement [in/in] =')
disp(max(abs(dtrn)))

figure('Name','RSN864090')
subplot(211); plot(tRSN864,VRSN864n)
xlabel('Time (sec)')
ylabel('Base Shear (k/k)')

subplot(212); plot(tRSN864,dtsn,tRSN864,dtrn)
xlabel('Time (sec)')
ylabel('Displacement (in/in)')
legend('SDOF','Roof')

figure('Name','RSN864090')
plot(maxdrift,hi/12)
xlabel('Story Drift [%]')
ylabel('Height [ft]')

%   RSN6915X1
run RSN6915_DARFIELD_HVSCS26W_a;
m = M1;
ag = column(ag);
p = m*ag*g;
c = 2*z*M1*w1;
k = M1*w1^2;
tRSN6915 = [0:dt:dt*(length(ag)-1)]';

% For Linear Response
gamma = 1/2;
beta = 1/4;
[aRSN6915,vRSN6915,uRSN6915] = NewmarkIntegrator(gamma,beta,m,c,k,p,dt);
aRSN6915 = abs(aRSN6915-p/m);

% Base Shear V(t) Normalized by W
VRSN6915n = M1*(w1^2)/W*uRSN6915;     % [k/k]

% Single DOF and Roof Level Displacements normalized by SDOF and Roof hts
d = zeros(length(hi),length(uRSN6915));

for i = 1:length(hi)
    d(i,:) = phi(i).*uRSN6915';
end

dtsn = uRSN6915 ./ (76.112*12);        % [in/in]
dtrn = d(9,:)' ./ (99*12);            % [in/in]

% Story Drift
drift = zeros([length(hi) length(uRSN6915)]);
maxdrift = zeros([length(hi) 1]);

for i = 2:length(hi)
    drift(i,:) = abs(d(i,:) - d (i-1,:))/(hi(i)-hi(i-1))*100; % [%]
    maxdrift(i) = max(drift(i,:));
end 

fprintf('\nPeak Response Quantities for RSN6915W are:\n')
disp('Normalized Base Shear [k/k] =')
disp(max(abs(VRSN6915n)))
disp('Base Shear [k] =')
disp(max(abs(VRSN6915n)*W))
disp('Overturning Moment [ft-k] =')
disp(max(abs(VRSN6915n)*W)*76.112)
disp('Normalized SDOF Displacement [in/in] =')
disp(max(abs(dtsn)))
disp('SDOF Displacement [in] =')
disp(max(abs(dtsn)*76.112*12))
disp('Normalized Roof Displacement [in/in] =')
disp(max(abs(dtrn)))

figure('Name','RSN6915W')
subplot(211); plot(tRSN6915,VRSN6915n)
xlabel('Time (sec)')
ylabel('Base Shear (k/k)')

subplot(212); plot(tRSN6915,dtsn,tRSN6915,dtrn)
xlabel('Time (sec)')
ylabel('Displacement (in/in)')
legend('SDOF','Roof')

figure('Name','RSN6915W')
plot(maxdrift,hi/12)
xlabel('Story Drift [%]')
ylabel('Height [ft]')

%   RSN1616X1
run RSN1616_DUZCE_362N_a;
m = M1;
ag = column(ag);
p = m*ag*g;
c = 2*z*M1*w1;
k = M1*w1^2;
tRSN1616 = [0:dt:dt*(length(ag)-1)]';

% For Linear Response
gamma = 1/2;
beta = 1/4;
[aRSN1616,vRSN1616,uRSN1616] = NewmarkIntegrator(gamma,beta,m,c,k,p,dt);
aRSN1616 = abs(aRSN1616-p/m);

% Base Shear V(t) Normalized by W
VRSN1616n = M1*(w1^2)/W*uRSN1616;     % [k/k]

% Single DOF and Roof Level Displacements normalized by SDOF and Roof hts
d = zeros(length(hi),length(uRSN1616));

for i = 1:length(hi)
    d(i,:) = phi(i).*uRSN1616';
end

dtsn = uRSN1616 ./ (76.112*12);        % [in/in]
dtrn = d(9,:)' ./ (99*12);            % [in/in]

% Story Drift
drift = zeros([length(hi) length(uRSN1616)]);
maxdrift = zeros([length(hi) 1]);

for i = 2:length(hi)
    drift(i,:) = abs(d(i,:) - d (i-1,:))/(hi(i)-hi(i-1))*100; % [%]
    maxdrift(i) = max(drift(i,:));
end 

fprintf('\nPeak Response Quantities for RSN1616N are:\n')
disp('Normalized Base Shear [k/k] =')
disp(max(abs(VRSN1616n)))
disp('Base Shear [k] =')
disp(max(abs(VRSN1616n)*W))
disp('Overturning Moment [ft-k] =')
disp(max(abs(VRSN1616n)*W)*76.112)
disp('Normalized SDOF Displacement [in/in] =')
disp(max(abs(dtsn)))
disp('SDOF Displacement [in] =')
disp(max(abs(dtsn)*76.112*12))
disp('Normalized Roof Displacement [in/in] =')
disp(max(abs(dtrn)))

figure('Name','RSN1616N')
subplot(211); plot(tRSN1616,VRSN1616n)
xlabel('Time (sec)')
ylabel('Base Shear (k/k)')

subplot(212); plot(tRSN1616,dtsn,tRSN1616,dtrn)
xlabel('Time (sec)')
ylabel('Displacement (in/in)')
legend('SDOF','Roof')

figure('Name','RSN1616N')
plot(maxdrift,hi/12)
xlabel('Story Drift [%]')
ylabel('Height [ft]')

%% Determine Peak Response Values based on 1st Mode Period and Spectral
%  ordinates from PEER Database

fprintf('\nBelow are the results of the Response Spectrum Analysis:\n')
RSN864AX1 = g*xlsread('_SearchResults.csv',1,'C236:C237'); % [in/sec2]
RSN1616AX1 = g*xlsread('_SearchResults.csv',1,'K236:K237'); % [in/sec2]
RSN6915AX1 = g*xlsread('_SearchResults.csv',1,'Z236:Z237');   % [in/sec2]

RSN864AX1 = (RSN864AX1(2)-RSN864AX1(1))/(1.1-1)*(T1-1)+RSN864AX1(1);
RSN1616AX1 = (RSN1616AX1(2)-RSN1616AX1(1))/(1.1-1)*(T1-1)+RSN1616AX1(1);
RSN6915AX1 = (RSN6915AX1(2)-RSN6915AX1(1))/(1.1-1)*(T1-1)+RSN6915AX1(1);

RSN864DX1 = (1/w1^2)*RSN864AX1;    % [in]
RSN1616DX1 = (1/w1^2)*RSN1616AX1;  % [in]
RSN6915DX1 = (1/w1^2)*RSN6915AX1;  % [in]

fprintf('\nPeak Response Quantities for RSN864090 are:\n')
disp('Base Shear [k] =')
disp(M1*(w1^2)*RSN864DX1)
disp('SDOF Displacement [in] =')
disp(RSN864DX1)
disp('Overturning Moment [ft-k] =')
disp(M1*(w1^2)*RSN864DX1*76.112)

fprintf('\nPeak Response Quantities for RSN6915W are:\n')
disp('Base Shear [k] =')
disp(M1*(w1^2)*RSN6915DX1)
disp('SDOF Displacement [in] =')
disp(RSN6915DX1)
disp('Overturning Moment [ft-k] =')
disp(M1*(w1^2)*RSN6915DX1*76.112)

fprintf('\nPeak Response Quantities for RSN1616N are:\n')
disp('Base Shear [k] =')
disp(M1*(w1^2)*RSN1616DX1)
disp('SDOF Displacement [in] =')
disp(RSN1616DX1)
disp('Overturning Moment [ft-k] =')
disp(M1*(w1^2)*RSN1616DX1*76.112)

%% Compute the Nonlinear Response of the System

z = 0.05;        % [% Damping]
M1 = 7.508;      % [k-sec^2/in]
W = (M1/0.67)*g; % [k]
hi = 12.*[0; 15; 27; 39; 51; 63; 75; 87; 99]; % [in]
phi = [0.000; 0.054; 0.164; 0.323; 0.518; 0.739; 0.977; 1.225; 1.476];
m = M1;
Fy = 0.15*W;
k = M1*w1^2;
ksh = 0.05*k;

%   RSN864090
run RSN864_LANDERS_JOS090_a;
ag = column(ag);
p = m*ag*g;
c = 2*z*M1*w1;
tRSN864 = [0:dt:dt*(length(ag)-1)]';

% For Non-Linear Response
gamma = 1/2;
beta = 1/4;
[aRSN864,vRSN864,uRSN864] = NewmarkIntegratorNL(gamma,beta,m,c,k,p,dt, ...
    Fy,ksh);
aRSN864 = abs(aRSN864-p/m);

% Base Shear V(t) Normalized by W
VRSN864n = M1*(w1^2)/W*uRSN864;     % [k/k]

% Single DOF and Roof Level Displacements normalized by SDOF and Roof hts
d = zeros(length(hi),length(uRSN864));

for i = 1:length(hi)
    d(i,:) = phi(i).*uRSN864';
end

dtsn = uRSN864 ./ (76.112*12);        % [in/in]
dtrn = d(9,:)' ./ (99*12);            % [in/in]

% Story Drift
drift = zeros([length(hi) length(uRSN864)]);
maxdrift = zeros([length(hi) 1]);

for i = 2:length(hi)
    drift(i,:) = abs(d(i,:) - d (i-1,:))/(hi(i)-hi(i-1))*100; % [%]
    maxdrift(i) = max(drift(i,:));
end 

fprintf('\nBelow are the results of the NL Response History Analysis:\n')
fprintf('\nPeak Response Quantities for RSN864090 are:\n')
disp('Normalized Base Shear [k/k] =')
disp(max(abs(VRSN864n)))
disp('Base Shear [k] =')
disp(max(abs(VRSN864n)*W))
disp('Overturning Moment [ft-k] =')
disp(max(abs(VRSN864n)*W)*76.112)
disp('Normalized SDOF Displacement [in/in] =')
disp(max(abs(dtsn)))
disp('SDOF Displacement [in] =')
disp(max(abs(dtsn)*76.112*12))
disp('Normalized Roof Displacement [in/in] =')
disp(max(abs(dtrn)))

figure('Name','RSN864090 - Nonlinear')
subplot(211); plot(tRSN864,VRSN864n)
xlabel('Time (sec)')
ylabel('Base Shear (k/k)')

subplot(212); plot(tRSN864,dtsn,tRSN864,dtrn)
xlabel('Time (sec)')
ylabel('Displacement (in/in)')
legend('SDOF','Roof')

figure('Name','RSN864090 - Nonlinear')
plot(maxdrift,hi/12)
xlabel('Story Drift [%]')
ylabel('Height [ft]')

%   RSN6915X1
run RSN6915_DARFIELD_HVSCS26W_a;
ag = column(ag);
p = m*ag*g;
c = 2*z*M1*w1;
tRSN6915 = [0:dt:dt*(length(ag)-1)]';

% For Non-Linear Response
gamma = 1/2;
beta = 1/4;
[aRSN6915,vRSN6915,uRSN6915] = NewmarkIntegratorNL(gamma,beta,m,c,k,p,...
    dt,Fy,ksh);
aRSN6915 = abs(aRSN6915-p/m);

% Base Shear V(t) Normalized by W
VRSN6915n = M1*(w1^2)/W*uRSN6915;     % [k/k]

% Single DOF and Roof Level Displacements normalized by SDOF and Roof hts
d = zeros(length(hi),length(uRSN6915));

for i = 1:length(hi)
    d(i,:) = phi(i).*uRSN6915';
end

dtsn = uRSN6915 ./ (76.112*12);        % [in/in]
dtrn = d(9,:)' ./ (99*12);            % [in/in]

% Story Drift
drift = zeros([length(hi) length(uRSN6915)]);
maxdrift = zeros([length(hi) 1]);

for i = 2:length(hi)
    drift(i,:) = abs(d(i,:) - d (i-1,:))/(hi(i)-hi(i-1))*100; % [%]
    maxdrift(i) = max(drift(i,:));
end 

fprintf('\nPeak Response Quantities for RSN6915W are:\n')
disp('Normalized Base Shear [k/k] =')
disp(max(abs(VRSN6915n)))
disp('Base Shear [k] =')
disp(max(abs(VRSN6915n)*W))
disp('Overturning Moment [ft-k] =')
disp(max(abs(VRSN6915n)*W)*76.112)
disp('Normalized SDOF Displacement [in/in] =')
disp(max(abs(dtsn)))
disp('SDOF Displacement [in] =')
disp(max(abs(dtsn)*76.112*12))
disp('Normalized Roof Displacement [in/in] =')
disp(max(abs(dtrn)))

figure('Name','RSN6915W - Nonlinear')
subplot(211); plot(tRSN6915,VRSN6915n)
xlabel('Time (sec)')
ylabel('Base Shear (k/k)')

subplot(212); plot(tRSN6915,dtsn,tRSN6915,dtrn)
xlabel('Time (sec)')
ylabel('Displacement (in/in)')
legend('SDOF','Roof')

figure('Name','RSN6915W - Nonlinear')
plot(maxdrift,hi/12)
xlabel('Story Drift [%]')
ylabel('Height [ft]')

%   RSN1616X1
run RSN1616_DUZCE_362N_a;
ag = column(ag);
p = m*ag*g;
c = 2*z*M1*w1;
tRSN1616 = [0:dt:dt*(length(ag)-1)]';

% For Non-Linear Response
gamma = 1/2;
beta = 1/4;
[aRSN1616,vRSN1616,uRSN1616] = NewmarkIntegratorNL(gamma,beta,m,c,k,p,...
    dt,Fy,ksh);
aRSN1616 = abs(aRSN1616-p/m);

% Base Shear V(t) Normalized by W
VRSN1616n = M1*(w1^2)/W*uRSN1616;     % [k/k]

% Single DOF and Roof Level Displacements normalized by SDOF and Roof hts
d = zeros(length(hi),length(uRSN1616));

for i = 1:length(hi)
    d(i,:) = phi(i).*uRSN1616';
end

dtsn = uRSN1616 ./ (76.112*12);        % [in/in]
dtrn = d(9,:)' ./ (99*12);            % [in/in]

% Story Drift
drift = zeros([length(hi) length(uRSN1616)]);
maxdrift = zeros([length(hi) 1]);

for i = 2:length(hi)
    drift(i,:) = abs(d(i,:) - d (i-1,:))/(hi(i)-hi(i-1))*100; % [%]
    maxdrift(i) = max(drift(i,:));
end 

fprintf('\nPeak Response Quantities for RSN1616N are:\n')
disp('Normalized Base Shear [k/k] =')
disp(max(abs(VRSN1616n)))
disp('Base Shear [k] =')
disp(max(abs(VRSN1616n)*W))
disp('Overturning Moment [ft-k] =')
disp(max(abs(VRSN1616n)*W)*76.112)
disp('Normalized SDOF Displacement [in/in] =')
disp(max(abs(dtsn)))
disp('SDOF Displacement [in] =')
disp(max(abs(dtsn)*76.112*12))
disp('Normalized Roof Displacement [in/in] =')
disp(max(abs(dtrn)))

figure('Name','RSN1616N - Nonlinear')
subplot(211); plot(tRSN1616,VRSN1616n)
xlabel('Time (sec)')
ylabel('Base Shear (k/k)')

subplot(212); plot(tRSN1616,dtsn,tRSN1616,dtrn)
xlabel('Time (sec)')
ylabel('Displacement (in/in)')
legend('SDOF','Roof')

figure('Name','RSN1616N - Nonlinear')
plot(maxdrift,hi/12)
xlabel('Story Drift [%]')
ylabel('Height [ft]')
