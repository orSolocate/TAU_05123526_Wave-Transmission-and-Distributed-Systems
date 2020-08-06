clc;
clear all;
close all;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%part B%%%%%%%%%%%%%%%%%%%%%%%
T = 1*10^-9; % sec
f0 = 10*10^9; % Hz
w0 = 2*pi*f0; % rad/s
n = 10^6; %number of samples
t = linspace(-2*T,2*T,n); % n time values
dt= (4*T)/(n-1); %time invervals 

F = heaviside(t+(T/2)) - heaviside(t-(T/2));%rectangular pulse
V_t = F.*sin(w0*t);

% time:
figure('Name','B1: VG(t)','NumberTitle','off');
p=plot(t/T,V_t,'-s');
p.MarkerSize=10;
p.MarkerIndices=3*n/8:(n/2)/20:5*n/8;
ylabel('Vg[V]');
xlabel('t/T []');
grid on;

%frequiency:
V_w=fft(V_t);
% converting the bins into [Hz] or [rad/sec]
k=1:n; %indices of f
Hz= (k-1)*(10^6*dt).^(-1); %equivalent Hz
w= Hz.*(2*pi);

figure('Name','B1: VG(w)','NumberTitle','off');
plot(w,abs(V_w));
axis([-10^11  10^11 0 10^5]) % zoom in on first data pair
ylabel('|F{Vg}|');
xlabel('w[red/sec]');
grid on;

%% Section C1
a = 20*10^-3;
c = 3*10^8; % speed of light m/s
f_cut=c/(2*a);
w_cut=2*pi*f_cut;
w_c=linspace(w_cut/c,3*w_cut/c,n);
k = sqrt((w_c).^2 -(pi/a)^2);%dispersion relation

figure('Name','C1: K(w) VS. w/c','NumberTitle','off');
plot(w_c,k);
xlabel("w/c");
ylabel('k(w/c)');

vp= c/sqrt(1-(f_cut/f0)^2)
vg= c*sqrt(1-(f_cut/f0)^2)

%% section C2
z_quarter = 1/(4*f0*(1/vg - 1/vp)) %using eq (9)

%% section C3
z = z_quarter.*[1,2,3,4,1/z_quarter,200/z_quarter];
F = @(t) heaviside(t+(T/2)) - heaviside(t-(T/2));
t=linspace(-20*T,20*T,n);
%used to compare with section D
figure('Name','C3: V(z,t) approximated','NumberTitle','off');

for i=1:length(z)
    tau=t-z(i)/vg;
    vzt_tilda=F(tau).*sin(w0*(tau+(z(i)/vg)-(z(i)/vp)));
    subplot(6,1,i);
    plot(tau,vzt_tilda);
    grid on;
    title(['V(z=',num2str(z(i)),',tau) approximated']);
    if z(i) ~= 200
        axis([-2*T 2*T -3 3]);
    end
    xlabel('tau [sec]');
    ylabel=('V(z,t) [V]');
end

%% section D
figure('Name','D: V(z,t) Exact','NumberTitle','off');
t=linspace(-20*T,20*T,n);
for i=1:length(z)
    tau=t-z(i)/vg;
    v_ztau= ifft(exp(1i.*w.*tau-1*1i*(k-w./vg).*z(i)).*V_w);
    subplot(6,1,i);
    plot(tau,v_ztau);
    grid on;
    title(['V(z=',num2str(z(i)),',tau) exact']);
    xlabel('tau [sec]');
    ylabel=('V(z,t) [v]');
end
