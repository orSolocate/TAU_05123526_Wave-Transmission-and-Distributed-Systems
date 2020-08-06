clc;
clear all;
close all;

%% Section B
%calculate central frequency
c=3*10^8; %[m/sec]
lambda1=400*10^-9; %[m]
lambda2=800*10^-9; %[m]
f1=c/lambda1;
f2=c/lambda2;
f0=(f1+f2)/2;

j=1i;
e0=8.85*10^-12; %[F/m]
er=2.25; %glass medium
eT=1.5; %transformer medium, from A
miu0=4*pi*10^-7; %[H/m]

eta_air=120*pi; %[ohm]
eta_glass=eta_air/sqrt(2.25); %[ohm]
eta_T=eta_air/sqrt(eT); %[ohm]
lambda_T0=c/(f0*sqrt(eT)); %[m]
dT=lambda_T0/4; %half wave tranformer [m]

f_vec=linspace(f1,f2,100);
sixty_in_rad = 60*pi/180;
theta_i_vec=linspace(-sixty_in_rad,sixty_in_rad,100);
[f,theta_i]=meshgrid(f_vec,theta_i_vec);

gemmaIn_TE_f_teta= GetGemmaIn(theta_i,f,'TE',eT,eta_T,dT,j,c,eta_air,eta_glass);
gemmaIn_TE_f_teta_0 = GetGemmaIn(0,f_vec,'TE',eT,eta_T,dT,j,c,eta_air,eta_glass);
gemmaIn_TE_f0_teta= GetGemmaIn(theta_i_vec,f0,'TE',eT,eta_T,dT,j,c,eta_air,eta_glass);

gemmaIn_TM_f_teta= GetGemmaIn(theta_i,f,'TM',eT,eta_T,dT,j,c,eta_air,eta_glass);
gemmaIn_TM_f_teta_0 = GetGemmaIn(0,f_vec,'TM',eT,eta_T,dT,j,c,eta_air,eta_glass);
gemmaIn_TM_f0_teta= GetGemmaIn(theta_i_vec,f0,'TM',eT,eta_T,dT,j,c,eta_air,eta_glass);

getPlots(gemmaIn_TE_f_teta,gemmaIn_TE_f_teta_0,gemmaIn_TE_f0_teta,gemmaIn_TM_f_teta,gemmaIn_TM_f_teta_0,gemmaIn_TM_f0_teta,theta_i_vec,theta_i,f,f_vec,1)

%% Section D
Eps={[1.257,1.773],[1.131,1.493,1.970],[1.0682,1.301,1.710,2.085]};
N=[2,3,4];
for i=1:length(Eps)
    eps = Eps{i};
    [d, etas]= GetLengths_and_etas(eps,eta_air);
    gemmaIn_TE_f_teta= GetGemmaIn(theta_i,f,'TE',eps,etas,d,j,c,eta_air,eta_glass);
    gemmaIn_TE_f_teta_0 = GetGemmaIn(0,f_vec,'TE',eps,etas,d,j,c,eta_air,eta_glass);
    gemmaIn_TE_f0_teta= GetGemmaIn(theta_i_vec,f0,'TE',eps,etas,d,j,c,eta_air,eta_glass);

    gemmaIn_TM_f_teta= GetGemmaIn(theta_i,f,'TM',eps,etas,d,j,c,eta_air,eta_glass);
    gemmaIn_TM_f_teta_0 = GetGemmaIn(0,f_vec,'TM',eps,etas,d,j,c,eta_air,eta_glass);
    gemmaIn_TM_f0_teta= GetGemmaIn(theta_i_vec,f0,'TM',eps,etas,d,j,c,eta_air,eta_glass);
    getPlots(gemmaIn_TE_f_teta,gemmaIn_TE_f_teta_0,gemmaIn_TE_f0_teta,gemmaIn_TM_f_teta,gemmaIn_TM_f_teta_0,gemmaIn_TM_f0_teta,theta_i_vec,theta_i,f,f_vec,N(i))
end
%% functions:
function [d,etas]= GetLengths_and_etas(eps,eta_air)
    c=3*10^8; %[m/sec]
    lambda1=400*10^-9; %[m]
    lambda2=800*10^-9; %[m]
    f1=c/lambda1;
    f2=c/lambda2;
    f0=(f1+f2)/2;
    lambda_T0=c./(f0.*sqrt(eps)); %[m]
    d=lambda_T0./4; %half wave tranformer [m]
    etas=eta_air./sqrt(eps);
end

function [gemmaIn]= GetGemmaIn(theta_i,f,polarization,EPS,ETAS,D,j,c,eta_air,eta_glass)
    n_glass =sqrt(2.25);
    theta_glass=asin(sin(theta_i)./n_glass); %snell law
    if polarization == 'TE'
        z_air=eta_air./cos(theta_i);
        z_glass=eta_glass./cos(theta_glass);
    else
        z_air=eta_air.*cos(theta_i);
        z_glass=eta_glass.*cos(theta_glass);
    end
    Zin=0;
    for i= length(ETAS):-1:1
        d= D(i);
        eta= ETAS(i); 
        theta=asin(sin(theta_i)./sqrt(EPS(i))); %snell law
        K= 2*pi.*f.*sqrt(EPS(i))./c; %[1/m]
        kz=K.*cos(theta); %[1/m]  
        if polarization == 'TE'
            Zc =eta./cos(theta);
        else
            Zc =eta.*cos(theta);
        end
        if i==length(ETAS)
            ZL=z_glass;
        end
        if i ==1
            Zin= Zc.*(ZL+j*Zc.*tan(kz.*d))./(Zc+j*ZL.*tan(kz.*d));
        else
            ZL= Zc.*(ZL+j*Zc.*tan(kz.*d))./(Zc+j*ZL.*tan(kz.*d));
        end
      end
    gemmaIn=(Zin-z_air)./(Zin+z_air);
end


function []= getPlots(gemmaIn_TE_f_teta,gemmaIn_TE_f_teta_0,gemmaIn_TE_f0_teta,gemmaIn_TM_f_teta,gemmaIn_TM_f_teta_0,gemmaIn_TM_f0_teta,theta_i_vec,theta_i,f,f_vec,N)
    %% TE PLOTS
    if N==1
        TITLE = sprintf('PLOTS FOR SECTION C');
    else
        TITLE = sprintf('PLOTS FOR %d LAYER', N);
    end
    
    figure('Position', [10 10 900 1600],'Name',TITLE,'NumberTitle','off');
    subplot(3,2,1);
    surf(theta_i,f,abs(gemmaIn_TE_f_teta));
    colorbar;
    title("|Gamma in TE mode|");
    xlabel('theta_i [rad]');
    ylabel('frequency [Hz]');

    subplot(3,2,3);
    plot(f_vec,abs(gemmaIn_TE_f_teta_0));
    title("|Gamma in TE mode|, theta=0");
    xlabel('frequency [Hz]');
    ylabel('|Gama in|');

    subplot(3,2,5);
    plot(theta_i_vec,abs(gemmaIn_TE_f0_teta));
    title("|Gamma in TE mode|, f=f0");
    xlabel('teta_i [rad]');
    ylabel('|Gama in|');

%% TM PLOTS
%     figure('Name',TM_NAME,'NumberTitle','off');
    subplot(3,2,2);
    surf(theta_i,f,abs(gemmaIn_TM_f_teta));
    colorbar;
    title("|Gamma in TM mode|");
    xlabel('theta_i [rad]');
    ylabel('frequency [Hz]');
    % zlim([0 1])

    subplot(3,2,4);
    plot(f_vec,abs(gemmaIn_TM_f_teta_0),'m');
    title("|Gamma in TM mode|, theta=0");
    xlabel('frequency [Hz]');
    ylabel('|Gama in|');

    subplot(3,2,6);
    plot(theta_i_vec,abs(gemmaIn_TE_f0_teta));
    hold on
    plot(theta_i_vec,abs(gemmaIn_TM_f0_teta),'m');
    title("|Gamma in TM mode|, f=f0");
    xlabel('teta_i [rad]');
    ylabel('|Gama in|');
    legend('TE','TM')
end