% CEE 246 Homework 3 - Nonlinear SDOF Response
%
% Ray Abbiatici
% Version 1.0/RJA/31-Jan-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

clear
clc

%% Select Max Response Motions from Chi Chi & Miyagi Data

T1 = 1.037;       % [sec]
w1 = (2*pi/T1); % [Hz]
g = 386.4;        % [in/sec2]

% Chi Chi - Note Horizontal 1 = E, Horizontal 2 = N
T = xlsread('ChiChi_SearchResults.csv',1,'A161:A271'); % [sec]
ChiChiEast = xlsread('ChiChi_SearchResults.csv',1,'B161:B271');
ChiChiNorth = xlsread('ChiChi_SearchResults.csv',1,'C161:C271');

% Miyagi - Note Horizontal 1 = NS, Horizontal 2 = EW
MiyagiNS = xlsread('Miyagi_SearchResults.csv',1,'B161:B271');
MiyagiEW = xlsread('Miyagi_SearchResults.csv',1,'C161:C271');

figure
plot(T,ChiChiEast,T,ChiChiNorth,T,MiyagiNS,T,MiyagiEW)
xlim([0 3])
xline(T1)
xlabel('T [sec]')
ylabel('pSa (g)')
legend('ChiChiEast','ChiChiNorth','MiyagiNS','MiyagiEW')

fprintf('The Controlling Ground Motions are ChiChiEast and MiyagiNS')

%% Plot Nonlinear (elastic, perfectly plastic) constant strength spectra...
%   for periods between 0.1 and 3.0 seconds, dt = 0.05 sec, Fy = CyW
%   for Cy = 0.1, 0.15, 0.25. Plot spectra along with mean and st. dev
%   for each strength ratio.

% System Properties
gamma = 1/2;
beta = 1/4;
Tn1 = 1.037;         % [sec]
z = 0.05;            % [%]
m = 7.508;           % [k-sec^2/in]
W = m*g/0.67;        % [k]
k = w1^2*m;          % [k-in]
h = 76.1112*12;      % [in]

% Spectral Range
T = [0.1:0.05:3];
T = T';

spectLan = zeros(length(T),1);
spectDar = zeros(length(T),1);
spectDuz = zeros(length(T),1);
spectChi = zeros(length(T),1);
spectMiy = zeros(length(T),1);
sd  = zeros(length(T),1);

%Strength Ratios
Cy = [0.1;0.15;0.25;];
table = zeros(length(Cy),7);

% Perform Newmark Nonlinear Integration
for i = 1:length(Cy)
    Fy = Cy(i)*W;
    for j = 1:length(T)
        
        w = 2*pi/(T(j));      % [Hz]
        k = w^2*m;            % [k/in]
        dy = Fy/k;            % [in]  
        c = 2*z*m*w;          % 
        ksh = 0;         % [k/in]

        run RSN864_LANDERS_JOS090_a;
        ag = column(ag);
        p_Landers090 = m*ag*g*1.5779;
        [aLan,~,~]=NewmarkIntegratorNL(gamma,beta,m,c,k,p_Landers090,dt,Fy,ksh);
        aLan = abs(aLan - p_Landers090/m);
        ddLan = aLan/w^2/dy;
        spectLan(j) = max(ddLan);
        
        run RSN6915_DARFIELD_HVSCS26W_a;
        ag = column(ag);
        p_DarfieldW = m*ag*g*1.4432;
        [aDar,~,~]=NewmarkIntegratorNL(gamma,beta,m,c,k,p_DarfieldW,dt,Fy,ksh);
        aDar = abs(aDar - p_DarfieldW/m);
        ddDar = aDar/w^2/dy;
        spectDar(j) = max(ddDar);

        run RSN1616_DUZCE_362N_a;
        ag = column(ag);
        p_DuzceN = m*ag*g*12.8769;
        [aDuz,~,~]=NewmarkIntegratorNL(gamma,beta,m,c,k,p_DuzceN,dt,Fy,ksh);
        aDuz = abs(aDuz - p_DuzceN/m);
        ddDuz = aDuz/w^2/dy;
        spectDuz(j) = max(ddDuz);

        run RSN1205_CHICHI_CHY041_E;
        ag = column(ag);
        p_ChiChiEast = m*ag*g*1.7205;
        [aChi,~,~]=NewmarkIntegratorNL(gamma,beta,m,c,k,p_ChiChiEast,dt,Fy,ksh);
        aChi = abs(aChi - p_ChiChiEast/m);
        ddChi = aChi/w^2/dy;
        spectChi(j) = max(ddChi);

        run RSN5776_IWATE_54010NS;
        ag = column(ag);
        p_MiyagiNS = m*ag*g*4.0271;
        [aMiy,v2,u2]=NewmarkIntegratorNL(gamma,beta,m,c,k,p_MiyagiNS,dt,Fy,ksh);
        aMiy = abs(aMiy - p_MiyagiNS/m);
        ddMiy = aMiy/w^2/dy;
        spectMiy(j) = max(ddMiy);
        
        all = [spectLan(j);spectDar(j);spectDuz(j);spectChi(j);...
            spectMiy(j)];

        sd(j) = std(all);

    end
    
    average = (spectLan+spectDar+spectDuz+spectChi+spectMiy)/5;
    sdplus = average+sd;
    sdminus = average-sd;

     w = 2*pi/(T1);      % [Hz]
     k = w1^2*m;         % [k/in]
     dy = Fy/k;          % [in]  

    table(i,1) = Cy(i);
    table(i,7) = (sdplus(20)-sdplus(19))/(1.05-1)*(T1-1)+sdplus(19);
    table(i,6) = (average(20)-average(19))/(1.05-1)*(T1-1)+average(19);
    table(i,5) = dy;
    table(i,4) = dy;
    table(i,3) = ((sdplus(20)-sdplus(19))/(1.05-1)*(T1-1)+sdplus(19))*dy;
    table(i,2) = ((average(20)-average(19))/(1.05-1)*(T1-1)+average(19))*dy;

    figure('Name',sprintf('Cy = %.2f',Cy(i)))
    plot (T,spectLan,T,spectDar,T,spectDuz,T,spectChi,T,spectMiy)
    hold on
    plot(T,average,T,sdplus,T,sdminus,'LineWidth',2);
    xlabel('T [sec]')
    ylabel('Displacement Ductility [in/in]')
    legend('Landers','Darfield','Duzce','Chi-Chi','Miyagi','Average',...
        '+1 St. Dev.','-1 St. Dev.' )

end

fprintf('\nThere is a notable decline in displacement ductility as period increases.')

%% Calculate constant strength spectra for strength ratios of Cy between
%   0.1 and 1.0 using an increment of 0.1. Store these values in a matrix
%   [U] with the strength ratio in the first column and the calculated
%   displacement ductility values in the remaining columns for each of
%   the following periods: T=0.1 to 3.0 for dT = 0.10 

Cy = [0.01:0.01:1.0]';
U = zeros(length(Cy),length(T)+1);

for i = 1:length(Cy)
    
    Fy = Cy(i)*W;
    U(i,1) = Cy(i);
        
    for j = 1:length(T)
        
            w = 2*pi/(T(j));      % [Hz]
            k = w^2*m;            % [k/in]
            dy = Fy/k;            % [in]  
            c = 2*z*m*w;          % 
            ksh = 0;              % [k/in]

            run RSN1205_CHICHI_CHY041_E;
            ag = column(ag);
            p_ChiChiEast = m*ag*g*1.7205;
            [aChi,~,~]=NewmarkIntegratorNL(gamma,beta,m,c,k,p_ChiChiEast,dt,Fy,ksh);
            aChi = abs(aChi - p_ChiChiEast/m);
            ddChi = aChi/w^2/dy;
            spectChi(j) = max(ddChi);
            U(i,j+1) =spectChi(j);
    end
end

e = 1;
f = 1;
t = 1;

for i = 1:length(Cy)
    for j = 2:length(T+1)
        if abs(U(i,j) - 8) < 1
            eight(e,1) = T(j-1);
            eight(e,2) = Cy(i);
            e = e+1;
        
        elseif abs(U(i,j) - 4) < 0.2
            four(f,1) = T(j-1);
            four(f,2) = Cy(i);
            f = f+1;
            
        elseif abs(U(i,j) - 2) < 0.1
            two(t,1) = T(j-1);
            two(t,2) = Cy(i);
            t = t+1;

        end
    end
end

eight = sortrows(eight);
four = sortrows(four);
two = sortrows(two);

figure('Name','Constant Ductility Spectra')
plot(two(:,1),two(:,2),four(:,1),four(:,2),eight(:,1),eight(:,2));
xlabel('T [sec]')
ylabel('Cy')
legend('2  in/in','4 in/in','8 in/in')

fprintf('\nLower displacement ductility values produce lower strength ratios.')

