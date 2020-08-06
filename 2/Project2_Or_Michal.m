clc;
clear all;
close all;
%% check_our_results
zL=2+1.6*1i;
d12= stub2(zL,1/8,'ss');

%% Section A
zC = 50;
Vp = 3*(10^8);  %m/s
f = 2*(10^9):10*(10^6):4*(10^9); %frequency range
f0=3*(10^9);
D = 0.0125;  %m
zL = 100+1j*60; %ohm
   
Lambda = Vp./f; %f %m
lambda0=Vp/f0; %central wavelength

beta =(2*pi)./Lambda;
beta0 =(2*pi)/lambda0;


%Option1
L_B1 = 0.4224*lambda0; %OUR ANSWER IS: 0.4224
L_B1 = 0.4539*lambda0; %STUB2 ANSWER IS: 0.4539

L_A1 = 0.4265*lambda0;  %OUR ANSWER IS: 0.4265
L_A1 = 0.425*lambda0;  %STUB2 ANSWER IS: 0.425

%Option2
L_B2 = 0.1045*lambda0;%OUR ANSWER IS: 0.1045
L_B2 = 0.101*lambda0; %STUB2 ANSWER IS: 0.101

L_A2 = 0.316*lambda0; 
L_A2 = 0.3269*lambda0;  %STUB2 ANSWER IS: 0.3269

%% Section B

GamaIn_1=simulate(zC,zL,beta,beta0,D, L_A1, L_B1);
GamaIn_2=simulate(zC,zL,beta,beta0,D,L_A2,L_B2);

figure('Name','Section B','NumberTitle','off');
subplot(2,1,2);
plot(f,abs(GamaIn_1));
grid on
ylabel('|Gamma_in|');
xlabel('Frequency [Hz]');
title('Double stub transmission matching - The second option = section B');

subplot(2,1,1);
plot(f,abs(GamaIn_2));
grid on
ylabel('|Gamma_in|');
xlabel('Frequency [Hz]');
title('Double stub transmission matching - The preferred option = section B');

%% Section C

GamaIn_1C=simulate(zC,zL,beta,beta0,D, L_A1+0.5*lambda0, L_B1);
GamaIn_2C=simulate(zC,zL,beta,beta0,D, L_A2+0.5*lambda0, L_B2);
 
figure('Name','Section C','NumberTitle','off');
subplot(2,1,2);
%plot(f,abs(GamaIn_1));
%hold on
plot(f,abs(GamaIn_1C));
grid on
%legend('LA','LA+0.5Lambda');
ylabel('|Gamma_in|');
xlabel('Frequency [Hz]');
title('section C - The second option');
subplot(2,1,1);
%plot(f,abs(GamaIn_2));
%hold on
plot(f,abs(GamaIn_2C));
grid on
%egend('LA','LA+0.5Lambda');
ylabel('|Gamma_in|');
xlabel('Frequency [Hz]');
title('section C - The preferred option');
 
%% Section D

[Pload1,Pin1]=calcPower(zC,GamaIn_1);
[Pload2, Pin2]=calcPower(zC,GamaIn_2);

figure('Name','Section D','NumberTitle','off');
subplot(2,1,2);
plot(f,Pload1);
grid on
ylabel('Pload');
xlabel('Frequency [Hz]');
title('Power absorbed by load - The second option');
subplot(2,1,1);
plot(f,Pload2);
grid on
ylabel('Pload');
xlabel('Frequency [Hz]');
title('Power absorbed by load - The preferred option');

 
%% Section D- explanation
    
figure('Name','Section D- exp','NumberTitle','off');
subplot(2,1,2);
plot(f,Pload1);
hold on
plot(f,Pin1);
hold on
plot(f,Pin1-Pload1);
legend('P load','P tot','P lost');
grid on
ylabel('Pload');
xlabel('Frequency [Hz]');
title('Power absorbed by load - The second option');
subplot(2,1,1);
plot(f,Pload2);
hold on
plot(f,Pin2);
hold on
plot(f,Pin2-Pload2);
legend('P load','P tot','P lost');
grid on
ylabel('Pload');
xlabel('Frequency [Hz]');
title('Power absorbed by load - The preferred option');

function [Pload,Ptot]=calcPower(zC,GamaIn)
    VG=1; %Volt
    ZG=75; %ohm
    
    Zin = zC.*(1+GamaIn)./(1-GamaIn);
    Vin = VG.*(Zin./(ZG+Zin)); % voltage divider
    Vin_plus = Vin./(1+GamaIn);
    P_plus = 0.5.*(1./zC).*(abs(Vin_plus).^2); % this is the power that "enters" the transmishion line
    Pload = P_plus.*(1- abs(GamaIn).^2);
    Ptot = P_plus;
end

function  GamaIn=simulate(zC,zL,beta,beta0,D,L_A,L_B)

    YL = 1./zL;
    YsA = 1./(1j*zC*tan(beta.*L_A));
    YA = YsA+YL;
    Z_A = 1./YA;
    GamaA = (Z_A-zC)./(Z_A+zC);

    %calculate YinB
    Z_B = zC*((1+GamaA.*exp(-2.*1j.*beta.*D))./(1-GamaA.*exp(-2*1j.*beta.*D)));
    YB = (1./Z_B); % Y3

    %calculate Yin
    YsB = 1./(1j*zC*tan(beta.*L_B));
    Yin = YsB + YB;
    Zin = 1./Yin;
    GamaIn = (Zin-zC)./(Zin+zC); % Gama In Calculation
end