%%
clear all;
close all;
clc;
%%
% Variables Setting for Part1
Zc=100; %[ohm]
v=3*(10^8); %[m/sec] 
L=1; %[m]
f=1*(10^9); %[Hz]
w=2*pi*f;
RG=200; %[ohm]
RL=50; %[ohm]
coef_G= (RG-Zc)/(RG+Zc);
coef_L=(RL-Zc)/(RL+Zc);
T=L/v; %[sec]
V_D= Zc/(RG+Zc); %voltage divder, relation between V(z,t) and V(t)
%%
% Vectors setting for section B
t=0:T/1000:(10*T); % Time vector
V_B=0;
z0=L/2;
for i = 1:5
   Vi_pls =V_D*(coef_L^(i-1))*(coef_G^(i-1))*sin(w*(t-(z0/v)-2*T*(i-1))).*heaviside((t-(z0/v)-2*T*(i-1))); % for the V+ vector
   Vi_min =V_D*(coef_L^(i))*(coef_G^(i-1))*sin(w*(t+(z0/v)-2*T*i)).*heaviside((t+(z0/v)-2*T*i)); % for the V- vector
   V_B=V_B+Vi_pls+Vi_min; % Summing the voltage vector
end

%% section (B)
figure('Name','part 1- B','NumberTitle','off');

subplot(3,1,1);
plot(t(1:2001),V_B(1:2001));
grid on
set(gca,'FontSize',12);
title ('V(L/2,t) , 0<=t<=2T');
xlabel('t');
ylabel('V(L/2,t)');

subplot(3,1,2);
plot(t(2001:4001),V_B(2001:4001));
grid on
set(gca,'FontSize',12);
title ('V(L/2,t) , 2T<=t<=4T');
xlabel('t');
ylabel('V(L/2,t)');

subplot(3,1,3);
plot(t(8001:10001),V_B(8001:10001));
grid on
set(gca,'FontSize',12);
title ('V(L/2,t) , 8T<=t<=10T');
xlabel('t');
ylabel('V(L/2,t)');


%% section (C)
z=0:L/1000:L; % space vector
times=[T/2,3*T/2,10*T];
V_C={};
for index=1:3
    V_C{index}=0;
    for i =1:(floor(times(index)/(2*T))+1)
       Vi_pls =V_D*(coef_L^(i-1))*(coef_G^(i-1))*sin(w*(times(index)-(z./v)-2*T*(i-1))).*heaviside((times(index)-(z/v)-2*T*(i-1))); % for the V+ vector
       Vi_min =V_D*(coef_L^(i))*(coef_G^(i-1))*sin(w*(times(index)+(z./v)-2*T*i)).*heaviside((times(index)+(z./v)-2*T*i)); % for the V- vector
       V_C{index}=V_C{index}+Vi_pls+Vi_min; % Summing the voltage vector
    end
end

figure('Name','part 1- C','NumberTitle','off');
subplot(3,1,1);
plot(z,V_C{1});
grid on
set(gca,'FontSize',12);
title ('V(z,t) , t=T/2');
xlabel('z');
ylabel('V(z,T/2)');

subplot(3,1,2);
plot(z,V_C{2});
grid on
set(gca,'FontSize',12);
title ('V(z,t) , t=3T/2');
xlabel('z');ylabel('V(z,3T/2)');

subplot(3,1,3);
plot(z,V_C{3});
grid on
set(gca,'FontSize',12);
title ('V(z,t) , t=10T');
xlabel('z');
ylabel('V(z,10T)');

%% section  (D)
beta= w/v;
Vg_p=-1i;
Zin_0=Zc*(1+coef_L*exp(-2*1i*beta*L))/(1-coef_L*exp(-2*1i*beta*L));
Vin_0_p=Vg_p*Zin_0/(Zin_0+RG);
V0_plus_p=Vin_0_p/(1+coef_L*exp(-2*1i*beta*L));
Vz=V0_plus_p.*exp(-1i*beta.*z).*(1+coef_L.*exp(-2*1i*beta.*(L-z))); % Vz calculation
Vz= transpose(Vz);
Vzt=real(Vz.*exp(1i.*w.*t)); % V_z.t calculation
%Vzt is a matrix; the rows rep' the z, the cols rep' the time

%% Part 2 - Section E-B
figure('Name','part 1-E.B','NumberTitle','off');

subplot(3,1,1);
plot(t(1:2001),Vzt(500,1:2001));
grid on
set(gca,'FontSize',12);
title ('V(L/2,t) , 0<=t<=2T');
xlabel('t');
ylabel('V(L/2,t)');

subplot(3,1,2);
plot(t(2001:4001),Vzt(500,2001:4001));
grid on
set(gca,'FontSize',12);
title ('V(L/2,t) , 2T<=t<=4T');
xlabel('t');
ylabel('V(L/2,t)');

subplot(3,1,3);
plot(t(8001:10001),Vzt(500,8001:10001));
grid on
set(gca,'FontSize',12);
title ('V(L/2,t) , 8T<=t<=10T');
xlabel('t');
ylabel('V(L/2,t)');
%% Part 2 - Section E-C

figure('Name','part 1-E.C','NumberTitle','off');
subplot(3,1,1);


plot(z,Vzt(:,500));
grid on
set(gca,'FontSize',12);
title ('V(z,t) , t=T/2');
xlabel('z');
ylabel('V(z,T/2)');

subplot(3,1,2);
plot(z,Vzt(:,1500));
grid on
set(gca,'FontSize',12);
title ('V(z,t) , t=3T/2');
xlabel('z');
ylabel('V(z,3T/2)');

subplot(3,1,3);
plot(z,Vzt(:,10001));
grid on
set(gca,'FontSize',12);
title ('V(z,t) , t=10T');
xlabel('z');
ylabel('V(z,10T)');

%% compare results

figure('Name','Compare results -Section B','NumberTitle','off');

subplot(3,1,1);
plot(t(1:2001),V_B(1:2001));
hold on;
plot(t(1:2001),Vzt(500,1:2001));
grid on;
set(gca,'FontSize',12);
title ('V(L/2,t) , 0<=t<=2T');
xlabel('t');
ylabel('V(L/2,t)');
legend('part1','part2');

subplot(3,1,2);
plot(t(2001:4001),V_B(2001:4001));
hold on;
plot(t(2001:4001),Vzt(500,2001:4001));
grid on;
set(gca,'FontSize',12);
title ('V(L/2,t) , 2T<=t<=4T');
xlabel('t');
ylabel('V(L/2,t)');
legend('part1','part2');

subplot(3,1,3);
plot(t(8001:10001),V_B(8001:10001));
hold on;
plot(t(8001:10001),Vzt(500,8001:10001));
grid on;
set(gca,'FontSize',12);
title ('V(L/2,t) , 8T<=t<=10T');
xlabel('t');
ylabel('V(L/2,t)');
legend('part1','part2');

figure('Name','Compare results -Section C','NumberTitle','off');

subplot(3,1,1);
plot(z,V_C{1});
hold on;
plot(z,Vzt(:,500));
grid on;
set(gca,'FontSize',12);
title ('V(z,t) , t=T/2');
xlabel('z');
ylabel('V(z,T/2)');
legend('part1','part2');

subplot(3,1,2);
plot(z,V_C{2});
hold on;
plot(z,Vzt(:,1500));
grid on;
set(gca,'FontSize',12);
title ('V(z,t) , t=3T/2');
xlabel('z');ylabel('V(z,3T/2)');
legend('part1','part2');

subplot(3,1,3);
plot(z,V_C{3});
hold on;
plot(z,Vzt(:,10001));
grid on;
set(gca,'FontSize',12);
title ('V(z,t) , t=10T');
xlabel('z');
ylabel('V(z,10T)');
legend('part1','part2');