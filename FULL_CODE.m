%% COMPLETE SCRIPT
%By: Jesus Ferrand

clc;
clear;
close all
%% ONX INPUTS
I = struct('A',0.3,... % Engine Bypass Ratio. [Dimensionless]
    'B',0.01,... % Customer Core Bleed Ratio. [Dimensionless]
    'm_dot',100,... % Engine Intake Mass Flow Rate. [lbm/s]
    'M_0',1.6,... % Freestream Mach Number. [Dimensionless]
    'M_6',0.4,... % Core Exit Mach Number. [Dimensionless]
    'day','Standard',... % Atmospheric condition setting. [N/A]
    'Alt',30000,... % Altitude. [ft]
    'Tt_4',3200,... % Turbine Exit stagnation Temperature (TIT). [R]
    'Tt_7',3600,... % Afterburner Exit stagnation Temperature. [R]
    'h_PR',18000,... % Fuel Heating Value. [Btu/lbm]
    'CTOL',0.01,... % Low-Pressure Power takeoff coefficient. [Dimensionless]
    'CTOH',0.00,... % High-Pressure Power takeoff coefficient. [Dimensionless]
    'e1',0.05,... % LPT core entrance cooling air fraction. [Dimensionless]
    'e2',0.05,... % HPT core entrance cooling air fraction. [Dimensionless]
    'e_f',0.890,... % Fan polytropic efficiency (e_f). [Dimensionless]
    'e_cL',0.890,... % LPC polytropic efficiency (e_cL). [Dimensionless]
    'e_cH',0.900,... % HPC polytropic efficiency (e_cH). [Dimensionless]
    'e_tH',0.890,... % HPT polytropic efficiency (e_tH). [Dimensionless]
    'e_tL',0.910,... % LPT polytropic efficiency (e_tL). [Dimensionless]
    'eta_b',0.98,... % Burner adiabatic efficiency (eta_b). [Dimensionless]
    'eta_AB',0.97,... % Afterburner adiabatic efficiency (eta_AB). [Dimensionless]
    'eta_mL',0.99,... % LP spool mechanical efficiency (eta_mL). [Dimensionless]
    'eta_mH',0.98,... % HP spool mechanical efficiency (eta_mH). [Dimensionless]
    'eta_mPH',0.98,... % HP PTO mechanical efficiency (eta_mPH). [Dimensionless]
    'eta_mPL',0.98,... % LP PTO mechanical efficiency (eta_mPL). [Dimensionless]
    'P0_P9',1,... % Overall Engine Expansion Ratio. [Dimensionless]
    'pi_AB',0.96,... % Afterburner total pressure recovery. [Dimensionless]
    'pi_f',3.5,... % Fan pressure ratio. [Dimensionless]
    'pi_c',16,... % Overall compressor pressure ratio (OPR). [Dimensionless]
    'pi_cL',3.5,... % LPC pressure ratio. [Dimensionless]
    'pi_b',0.97,... % Main Burner total pressure recovery. [Dimensionless]
    'pi_d_max',0.97,... % Maximum diffuser pressure ratio (pi_d_max). [Dimensionless]
    'pi_M_max',0.97,... % Maximum mixer pressure ratio (pi_M_max). [Dimensionless]
    'pi_n',0.98); % Nozzle total pressure recovery. [Dimensionless]
%% MISSION "LEGS" (MANEUVERS) & WEIGHT
%n: load factor
%N: Number of Turns
%i: steps
L = {struct('ID',9,'Name','Warm-Up',...
    'AB',0,'Alt_1',0,'Alt_2',0,'CDR',0,'CL_M',0,'d',0,'i',0,'k_TO',0,...
    'mu_TO',0,'M_1',0,'M_2',0,'n',1,'N',0,'t',60,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);...
    struct('ID',4,'Name','Takeoff Acceleration',...
    'AB',1,'Alt_1',0,'Alt_2',0,'CDR',0.07,'CL_M',2,'d',0,'i',0,'k_TO',1.2,...
    'mu_TO',0.05,'M_1',0.09,'M_2',0,'n',1,'N',0,'t',60,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);...
    struct('ID',10,'Name','Takeoff Rotation',...
    'AB',1,'Alt_1',0,'Alt_2',0,'CDR',0.05,'CL_M',2,'d',0,'i',0,'k_TO',1.2,...
    'mu_TO',0,'M_1',0,'M_2',0,'n',1,'N',0,'t',0,'t_R',3,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',2,'Name','Horizontal Acceleration',...
    'AB',1,'Alt_1',0,'Alt_2',0,'CDR',0,'CL_M',0,'d',0,'i',1,'k_TO',1.2,...
    'mu_TO',0,'M_1',0.181,'M_2',0.7,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',3,'Name','Climb & Acceleration',...
    'AB',1,'Alt_1',0,'Alt_2',16000,'CDR',0,'CL_M',0,'d',0,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',0.7,'M_2',0.85,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',3,'Name','Climb & Acceleration',...
    'AB',1,'Alt_1',16000,'Alt_2',30000,'CDR',0,'CL_M',0,'d',0,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',0.85,'M_2',0.9,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',3,'Name','Climb & Acceleration',...
    'AB',1,'Alt_1',30000,'Alt_2',40000,'CDR',0,'CL_M',0,'d',0,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',0.9,'M_2',0.9,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);    
    struct('ID',7,'Name','BCA/BCM',...
    'AB',0,'Alt_1',40000,'Alt_2',0,'CDR',0,'CL_M',0,'d',300*6080,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',0.9,'M_2',0.9,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',3,'Name','Climb & Acceleration',...
    'AB',0,'Alt_1',41000,'Alt_2',30000,'CDR',0,'CL_M',0,'d',0,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',0.9,'M_2',1.6,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',5,'Name','Constant Speed Cruise',...
    'AB',1,'Alt_1',30000,'Alt_2',0,'CDR',0,'CL_M',0,'d',100*6080,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',1.6,'M_2',0,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',12,'Name','Payload Delivery','WP',2343);
    struct('ID',6,'Name','Constant Altitude Turn',...
    'AB',1,'Alt_1',30000,'Alt_2',0,'CDR',0,'CL_M',0,'d',100*6080,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',1.6,'M_2',0,'n',5,'N',2,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',6,'Name','Constant Altitude Turn',...
    'AB',1,'Alt_1',30000,'Alt_2',0,'CDR',0,'CL_M',0,'d',100*6080,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',0.9,'M_2',0,'n',5,'N',4,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',2,'Name','Horizontal Acceleration',...
    'AB',1,'Alt_1',30000,'Alt_2',0,'CDR',0,'CL_M',0,'d',0,'i',1,'k_TO',1.2,...
    'mu_TO',0,'M_1',0.9,'M_2',1.6,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',5,'Name','Constant Speed Cruise',...
    'AB',1,'Alt_1',30000,'Alt_2',0,'CDR',0,'CL_M',0,'d',100*6080,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',1.6,'M_2',0,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',11,'Name','Constant Energy Height',...
    'AB',0,'Alt_1',30000,'Alt_2',49000,'CDR',0,'CL_M',0,'d',100*6080,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',1.6,'M_2',0.9,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'frac',0.5,'par2',0);
    struct('ID',7,'Name','BCA/BCM',...
    'AB',0,'Alt_1',49000,'Alt_2',0,'CDR',0,'CL_M',0,'d',300*6080,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',0.9,'M_2',0.9,'n',1,'N',0,'t',0,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);
    struct('ID',8,'Name','Loiter',...
    'AB',0,'Alt_1',20000,'Alt_2',0,'CDR',0,'CL_M',0,'d',300*6080,'i',1,'k_TO',0,...
    'mu_TO',0,'M_1',0.9,'M_2',0.9,'n',1,'N',0,'t',1200,'t_R',0,'V_1',0,...
    'V_2',0,'par1',0,'par2',0);};

Engine = 'LBTF';
TL = 1.25; % Aircraft Thrust Loading. [Dimensionless]
WL = 64; % Aircraft Wing Loading. [psf]
S = 375; % Aircraft Planform. [ft^2]
WE = 22166; % Aircraft empty weight. [lbf]
WTO = 42000;
%% Function Callout
%R = ONX_LBTF_VSH(I);
tic
MISSION(L,Engine,I,TL,WL,WTO)
toc
%% MISSION Function
%TO-DO: "Code" par1 and par 2 for landing and constant Ze
function MISSION(L,Engine,I,TL,WL,WTO)
f = figure;
ax1 = axes('Parent',f,'Box','on');
hold(ax1,'on')
ylabel(ax1,'Instantaneous Weight Fraction ($\beta = \frac{W}{W_{TO}}$)','interpreter','latex')
xlabel(ax1,'Mission Stage','interpreter','latex')
g_c = 32.174;
R = ONX_LBTF_VSH(I);
O_MIL = OFX_LBTF_VSH(R,I,0,0,0);
O_MAX = OFX_LBTF_VSH(R,I,0,0,1);
[a,~] = size(L);
data = zeros(a,4);
beta = 1;
for i = 1:a
    switch L{i}.ID
        case 1  %Constant Speed Climb (Type A)
            Wf_Wi = 1;
            beta_temp = beta;
            s = 0;
            t = 0;
            h = linspace(L{i}.Alt_1,L{i}.Alt_2,L{i}.N+1); %Height "Stations." [ft]
            for k = 1:L{i}.i
                [~,~,P_m,~] = ATMOSphere(h_m,'BE',I.day);
                [~,a_1,~,~] = ATMOSphere(h(k),'BE',I.day);
                [~,a_2,~,~] = ATMOSphere(h(k+1),'BE',I.day);
                M_1 = a_1*L{i}.V_1;
                M_2 = a_2*L{i}.V_1;
                M_m = (M_1+M_2)/2;
                [alpha,O] = POWER_TOGGLE(R,I,M_m,L{i}.h,L{i}.AB,O_MIL,O_MAX);
                DZE = h(k+1) - h(k);
                CL =2*beta*WL/(1.4*P_m*M_m^2);
                [K1,~,~,CD0] = K1_CD0(M_m);
                CD = K1*CL^2 + CD0 + L{i}.CDR; %Aircraft Coefficient of Drag. [Dimensionless]
                u = CD*beta/(TL*CL*alpha);
                if DZE < 0
                    Wf_Wi_temp = 1;
                else
                    ds = DZE*u*CL/((1-u)*CD);
                    dt = ds/L{i}.V_1;
                    s = s+ds;
                    t = t+dt;
                    Wf_Wi_temp = exp(-O.S*DZE/(V_m*(1-u)*3600)); %Weight Fraction across leg. [Dimensionless]
                end      
                Wf_Wi = Wf_Wi*Wf_Wi_temp;
                beta = beta*Wf_Wi_temp;
            end
            beta = beta_temp;
        case 2  %Horizontal Acceleration (Type A)
            Wf_Wi = 1;
            beta_temp = beta;
            s = 0;
            t = 0;
            M = linspace(L{i}.M_1,L{i}.M_2,L{i}.i+1); %Mach Number "Stations." [Dimensionless]
            for k = 1:L{i}.i
                M_m = (M(k)+M(k+1))/2;
                [alpha,O] = POWER_TOGGLE(R,I,M_m,L{i}.Alt_1,L{i}.AB,O_MIL,O_MAX);
                [~,a_m,P_m,~] = ATMOSphere(L{i}.Alt_1,'BE',I.day);
                Vi = M(k)*a_m;
                Vf = M(k+1)*a_m;
                DZE = (Vf^2-Vi^2)/(2*g_c);
                CL =2*beta*WL/(1.4*P_m*M_m^2);
                [K1,~,~,CD0] = K1_CD0(M_m);
                CD = K1*CL^2 + CD0 + L{i}.CDR; %Aircraft Coefficient of Drag. [Dimensionless]
                u = CD*beta/(TL*CL*alpha);
                V_m = M_m*a_m;
                if DZE < 0
                    Wf_Wi_temp = 1;
                else
                    ds = DZE*u*CL/((1-u)*CD);
                    dt = ds/V_m;
                    s = s+ds;
                    t = t+dt;
                    Wf_Wi_temp = exp(-O.S*DZE/(V_m*(1-u)*3600)); %Weight Fraction across leg. [Dimensionless]
                end      
                Wf_Wi = Wf_Wi*Wf_Wi_temp;
                beta = beta*Wf_Wi_temp;
            end
            beta = beta_temp;
        case 3  %Climb and Acceleration (Type A)
            Wf_Wi = 1;
            beta_temp = beta;
            s = 0;
            t = 0;
            M = linspace(L{i}.M_1,L{i}.M_2,L{i}.i+1); %Mach Number "Stations." [Dimensionless]
            h = linspace(L{i}.Alt_1,L{i}.Alt_2,L{i}.i+1); %Height "Stations." [ft]
            for k = 1:L{i}.i
                h_m = (h(k)+h(k+1))/2;
                M_m = (M(k)+M(k+1))/2;
                [alpha,O] = POWER_TOGGLE(R,I,M_m,h_m,L{i}.AB,O_MIL,O_MAX);
                [~,a_m,P_m,~] = ATMOSphere(h_m,'BE',I.day);
                [~,a_1,~,~] = ATMOSphere(h(k),'BE',I.day);
                [~,a_2,~,~] = ATMOSphere(h(k+1),'BE',I.day);
                Vi = M(k)*a_1;
                Vf = M(k+1)*a_2;
                DZE = h(k+1)-h(k) + (Vf^2-Vi^2)/(2*g_c);
                CL =2*beta*WL/(1.4*P_m*M_m^2);
                [K1,~,~,CD0] = K1_CD0(M_m);
                CD = K1*CL^2 + CD0 + L{i}.CDR; %Aircraft Coefficient of Drag. [Dimensionless]
                u = CD*beta/(TL*CL*alpha);
                V_m = M_m*a_m;
                if DZE < 0
                    Wf_Wi_temp = 1;
                else
                    ds = DZE*u*CL/((1-u)*CD);
                    dt = ds/V_m;
                    s = s+ds;
                    t = t+dt;
                    Wf_Wi_temp = exp(-O.S*DZE/(V_m*(1-u)*3600)); %Weight Fraction across leg. [Dimensionless]
                end      
                Wf_Wi = Wf_Wi*Wf_Wi_temp;
                beta = beta*Wf_Wi_temp;
            end
            beta = beta_temp;
        case 4  %Takeoff Acceleration (Type A)
            [~,a_0,~,rho_0] = ATMOSphere(L{i}.Alt_1,'BE',I.day);
            VTO = L{i}.k_TO*sqrt(2*beta*WL/(rho_0*L{i}.CL_M)); %Aircraft Takeoff Speed (1.2*V_stall). [ft/s]
            M_TO = VTO/a_0;
            [alpha,O] = POWER_TOGGLE(R,I,M_TO,L{i}.Alt_1,L{i}.AB,O_MIL,O_MAX);
            [K1,~,~,CD0] = K1_CD0(M_TO);
            E_TO = L{i}.CDR + CD0 + K1*(L{i}.CL_M/(L{i}.k_TO^2))^2 - L{i}.mu_TO*L{i}.CL_M/(L{i}.k_TO^2);
            q = 0.5*rho_0*VTO^2; %Dynamic Pressure. [psf]
            u = (E_TO*q/(beta*WL) + L{i}.mu_TO)*(beta/(alpha*TL));
            Wf_Wi = exp(-O.S*VTO/(g_c*(1-u)*3600)); %Weight Fraction across leg. [Dimensionless]
        case 5  %Constant Altitude/Speed Cruise (Type B)
            O = OFX_LBTF_VSH(R,I,L{i}.M_1,L{i}.Alt_1,L{i}.AB);
            L{i}.V_1 = L{i}.M_1*O.a_0;
            q = 0.5*O.rho_0*(L{i}.V_1^2);
            [K1,~,~,CD0] = K1_CD0(L{i}.M_1);
            CL = beta*WL/q;
            CD = K1*CL^2+CD0;
            Wf_Wi = exp(-O.S*(CD+L{i}.CDR)*L{i}.d/(L{i}.V_1*CL*3600));
        case 6  %Constant Altitude/Speed Turn (Type B)
            O = OFX_LBTF_VSH(R,I,L{i}.M_1,L{i}.Alt_1,L{i}.AB);
            L{i}.V_1 = L{i}.M_1*O.a_0;
            q = 0.5*O.rho_0*(L{i}.V_1^2);
            CL = L{i}.n*beta*WL/q;
            [K1,~,~,CD0] = K1_CD0(L{i}.M_1);
            CD = K1*CL^2 + CD0;
            Wf_Wi = exp(-O.S*(CD + L{i}.CDR)*L{i}.n*2*pi*L{i}.N*L{i}.V_1/(CL*g_c*sqrt(L{i}.n^2-1)*3600)); %Weight Fraction across leg. [Dimensionless]
        case 7  %Best Subsonic Cruise Mach Number and Altitude (BCM/BCA) (Type B)
            O = OFX_LBTF_VSH(R,I,L{i}.M_1,L{i}.Alt_1,L{i}.AB);
            [K1,~,~,CD0] = K1_CD0(L{i}.M_1); 
            Wf_Wi = exp(-O.S*L{i}.d*sqrt(4*(CD0+L{i}.CDR)*K1)/(L{i}.M_1*O.a_0*3600));
        case 8  %Subsonic Loiter (Type B)
            O = OFX_LBTF_VSH(R,I,L{i}.M_1,L{i}.Alt_1,L{i}.AB);
            [K1,~,~,CD0] = K1_CD0(L{i}.M_1);
            q = 0.5*O.rho_0*O.V_0^2;
            CL = beta*WL/q;
            CD = K1*CL^2+CD0+L{i}.CDR;
            Wf_Wi = exp(-O.S*L{i}.t*CD/(CL*3600));
        case 9  %Warm-Up (Type B)
            [alpha,O] = POWER_TOGGLE(R,I,L{i}.M_1,L{i}.Alt_1,L{i}.AB,O_MIL,O_MAX);
            Wf_Wi = 1 - O.S*TL*L{i}.t*alpha/(beta*3600);
        case 10 %Takeoff Rotation (Type B)
            [alpha,O] = POWER_TOGGLE(R,I,M_TO,L{i}.Alt_1,L{i}.AB,O_MIL,O_MAX);
            Wf_Wi = 1 - O.S*alpha*TL*L{i}.t_R/(beta*3600); %Weight Fraction across leg. [Dimensionless]
        case 11 %Constant Energy Height Maneuver (Type B)                
            h_m = (L{i}.Alt_1 + L{i}.Alt_2)/2;
            [~,a_m,~,rho_m] = ATMOSphere(h_m,'BE',I.day);
            [~,a_1,~,~] = ATMOSphere(L{i}.Alt_1,'BE',I.day);
            Ze_1 = L{i}.Alt_1 + ((L{i}.M_1*a_1)^2)/(2*g_c);
            V_m = sqrt(2*g_c*(Ze_1-h_m));
            M_m = V_m/a_m;
            O = OFX_LBTF_VSH(R,I,M_m,h_m,L{i}.AB);
            q = 0.5*rho_m*V_m^2;
            CL = beta*WL/q;
            [K1,~,~,CD0] = K1_CD0(L{i}.M_1);
            CD = K1*CL^2+CD0+L{i}.CDR;
            L{i}.t = abs(L{i}.Alt_2-L{i}.Alt_1)/(V_m*L{i}.frac);
            Wf_Wi = exp(-O.S*L{i}.t*CD/(CL*3600));
        case 12 %Payload
            Wf_Wi = (WTO - L{i}.WP)/WTO;
    end
    beta = beta*Wf_Wi; %Updated Instantaneous weight fraction. [Dimensionless]
    fuel_burnt = WTO*(1-Wf_Wi);
    WTO = WTO - fuel_burnt;
    fprintf('Leg %-2i (%-25s):Wf_Wi = %7.5f beta = %7.5f fuel burnt = %7.5f\n',i,L{i}.Name,Wf_Wi,beta,fuel_burnt) %Command Line Output of weight fractions.
    data(i,1) = i; %Store leg ID.
    data(i,2) = beta; %Store Weight Fraction.
    data(i,3) = Wf_Wi; %Store Maneuver Weight Ratio.
    data(i,4) = fuel_burnt; %Maneuver Fuel Burn.
    plot(ax1,[i-0.5,i-0.5],[0,beta],'Color',[0 0.5 1],'Linewidth',15)
end
    function [alpha,O] = POWER_TOGGLE(R,I,M_0,Alt,AB,O_MIL,O_MAX)
        O = OFX_LBTF_VSH(R,I,M_0,Alt,AB);
        if AB == 0
            alpha = O.F/O_MIL.F;
        else
            alpha = O.F/O_MAX.F;
        end
    end
end
%% ONX Calculator
%TO-DO
%- Collect all thermo outputs.

%--------------------------------------------------------------------
%Title: "Mattingly's LBTF Parametric Cycle Analysis w/ VSH in MATLAB"
%By: Jesus Ferrand
%Credits: Jack D. Mattingly
%--------------------------------------------------------------------
%INFO:
function R = ONX_LBTF_VSH(I)
R = struct('A',I.A,'A_prime',0,'A4_A45',0,'A45_A6',0,'A16_A6',0,'A8_A6_dry',0,...
    'A6_A6A',0,'A9_A8_off',0,'A9_A8_on',0,...
    'eta_cH',0,'eta_cL',0,'eta_tH',0,'eta_tL',0,'eta_f',0,'eta_P_on',0,...
    'eta_P_off',0,'eta_TH_on',0,'eta_TH_off',0,...
    'f',0,'f_AB',0,'f_6A',0,'f_0_on',0,'f_0_off',0,'F_m_dot_on',0,...
    'F_m_dot_off',0,'MFP_4',0,'MFP_45',0,'MFP_6',0,'MFP_6A',0,...
    'pi_ABdry',0,'pi_AB',I.pi_AB,'PtP_9_off',0,'PtP_9_on',0,...
    'pi_c',I.pi_c,'pi_cL',I.pi_cL,'pi_cH',0,'pi_d',0,'pi_f',I.pi_f,...
    'pi_M',0,'pi_M_ideal',0,'pi_r',0,'pi_tH',0,'pi_tL',0,...
    'm_dot',I.m_dot,'M_4',1,'M_45',1,'M_6',I.M_6,'M_16',0,'M_6A',0,...
    'M_8_off',0,'M_8_on',0,'M_9_on',0,'M_9_off',0,...
    'P_0',0,...
    'h_0',0,...
    'Pt16_Pt6',0,...
    'S_on',0,...
    'Cp_6A',0,'gamma_6A',0,...
    'S_off',0,'tau_cH',0,'tau_cL',0,'tau_f',0,'tau_m1',0,'tau_m2',0,...
    'tau_M',0,'tau_tH',0,'tau_tL',0,'tau_lambda',0,'tau_r',0,'theta_0',0,...
    'T9_T0_off',0,'T9_T0_on',0,'V9_V0_off',0,'V9_V0_on',0,'M9_M0_off',0,'M9_M0_on',0,...
    'P_TOL',0,'P_TOH',0,'T_0',0,'Tt_0',0,'Tt_25',0,'Tt_3',0,'Tt_4',I.Tt_4,...
    'Tt_45',0,'Tt_5',0,'Tt_6',0,'Tt_6A',0,'Tt_7',I.Tt_7,'Tt_13',0,'Tt_16',0,...
    'e1',I.e1,'e2',I.e2,'B',I.B,...
    'lmbme1me2',1-I.B-I.e1-I.e2);
g_c = 32.174;
R.pi_ABdry = 1-(1-I.pi_AB)/2;
[R.T_0,~,R.P_0,rho_0] = ATMOSphere(I.Alt,'BE',I.day);
[~,R.h_0,Pr_0,Phi_0,Cp_0,R_0,gamma_0,a_0] = FAIR(1,R.T_0,0,'BE');
R.P_TOL = I.CTOL*I.m_dot*R.h_0; %Low-Pressure Power Off-Take [Btu/s]
R.P_TOH = I.CTOH*I.m_dot*R.h_0; %High-Pressure Power Off-Take [Btu/s]
V_0 = I.M_0*a_0; %ONX Ambient Velocity. [ft/s]
ht_0 = R.h_0 + (V_0^2)/(2*g_c*778.16); %Stagnation Enthalpy at Ambient.
R.tau_r = ht_0/R.h_0; %Ram Heating.
R.theta_0 = R.tau_r*R.T_0/518.67; %ONX Throttle Ratio.
[R.Tt_0,ht_0,Prt_0,Phit_0,~,~,~,~] = FAIR(2,ht_0,0,'BE');
R.pi_r = Prt_0/Pr_0; %Ram Pressure Ratio .
eta_R_spec = MIL_E_5008B(I.M_0);
R.pi_d = I.pi_d_max*eta_R_spec; %Diffuser Pressure Ratio.
ht_2 = ht_0; %Stagnation Enthalpy entering Fan. (Adiabatic Assumption)
Prt_2 = Prt_0; %Reduced Pressure entering Fan. (Adiabatic Assumption)
Prt_13 = Prt_2*power(I.pi_f,1/I.e_f);
[R.Tt_13,ht_13,~,Phit_13,~,~,~,~] = FAIR(3,Prt_13,0,'BE');
R.tau_f = ht_13/ht_2; %Fan Enthalpy ratio.
Prt_13i = Prt_2*I.pi_f;
[~,ht_13i,~,~,~,~,~,~] = FAIR(3,Prt_13i,0,'BE');
R.eta_f = (ht_13i-ht_2)/(ht_13-ht_2);
Prt_25 = Prt_2*power(I.pi_cL,1/I.e_cL);
[R.Tt_25,ht_25,~,Phit_25,~,~,~,~] = FAIR(3,Prt_25,0,'BE');
R.tau_cL = ht_25/ht_2;
Prt_25i = Prt_2*I.pi_cL;
[~,ht_25i,~,~,~,~,~,~] = FAIR(3,Prt_25i,0,'BE');
R.eta_cL = (ht_25i-ht_2)/(ht_25-ht_2);
R.pi_cH = I.pi_c/I.pi_cL;
Prt_3 = Prt_25*power(R.pi_cH,1/I.e_cH);
[R.Tt_3,ht_3,~,Phit_3,~,~,~,~] = FAIR(3,Prt_3,0,'BE');
R.tau_cH = ht_3/ht_25;
Prt_3i = Prt_25*R.pi_cH;
[~,ht_3i,~,~,~,~,~,~] = FAIR(3,Prt_3i,0,'BE');
R.eta_cH = (ht_3i-ht_25)/(ht_3-ht_25);
R.f = FUEL_NR(ht_3,I.eta_b,I.h_PR,I.Tt_4);
lmbmelme2tlpf = (1-R.B-R.e1-R.e2)*(1+R.f); %Optimization Step.
[~,ht_4,Prt_4,Phit_4,~,~,~,~] = FAIR(1,R.Tt_4,R.f,'BE');
[TtT4,PtP4,R.MFP_4] = MASSFP_NR(R.Tt_4,R.f,R.M_4);
R.tau_lambda = ht_4/R.h_0;
R.tau_m1 = (lmbmelme2tlpf+I.e1*R.tau_r*R.tau_cL*R.tau_cH/R.tau_lambda)/(lmbmelme2tlpf+I.e1);
R.tau_tH = 1 - (R.tau_r*R.tau_cL*(R.tau_cH-1)+(1+I.A)*I.CTOH/I.eta_mPH)/((I.eta_mH*R.tau_lambda)*lmbmelme2tlpf+I.e1*R.tau_r*R.tau_cL*R.tau_cH/R.tau_lambda);
ht_41 = ht_4*R.tau_m1;
f_41 = R.f/(1+R.f+I.e1/R.lmbme1me2);
[Tt_41,~,Prt_41,Phit_41,~,~,~,~] = FAIR(2,ht_41,f_41,'BE');
ht_44 = ht_41*R.tau_tH;
[Tt_44,~,Prt_44,Phit_44,~,~,~,~] = FAIR(2,ht_44,f_41,'BE');
R.pi_tH = power(Prt_44/Prt_41,1/I.e_tH);
Prt_44i = R.pi_tH*Prt_41;
[~,ht_44i,~,~,~,~,~,~] = FAIR(3,Prt_44i,f_41,'BE');
R.eta_tH = (ht_41-ht_44)/(ht_41-ht_44i);
R.tau_m2 = (lmbmelme2tlpf+I.e1+I.e2*(R.tau_r*R.tau_cL*R.tau_cH/(R.tau_lambda*R.tau_m1*R.tau_tH)))/(lmbmelme2tlpf+I.e1+I.e2);
ht_45 = ht_44*R.tau_m2;
f_45 = R.f/(1+R.f+(I.e1+I.e2)/R.lmbme1me2);

[R.Tt_45,~,Prt_45,Phit_45,~,~,~,~] = FAIR(2,ht_45,f_45,'BE');
[TtT45,PtP45,R.MFP_45] = MASSFP_NR(R.Tt_45,f_45,R.M_45);
R.A4_A45 = R.pi_tH*sqrt(R.Tt_4/R.Tt_45)*(R.MFP_45/R.MFP_4)/(1+(I.e1+I.e2)/lmbmelme2tlpf);
R.tau_tL = 1 - (R.tau_r*((R.tau_cL-1)+I.A*(R.tau_f-1))+(1+I.A)*I.CTOL/I.eta_mPL)/  (I.eta_mL*R.tau_lambda*R.tau_tH*(lmbmelme2tlpf+(I.e1+I.e2/R.tau_tH)*(R.tau_r*R.tau_cL*R.tau_cH/R.tau_lambda)));
ht_5 = ht_45*R.tau_tL;
[R.Tt_5,~,Prt_5,Phit_5,~,~,~,~] = FAIR(2,ht_5,f_45,'BE');
R.pi_tL = power(Prt_5/Prt_45,1/I.e_tL);
Prt_5i = R.pi_tL*Prt_45;
[~,ht_5i,~,~,~,~,~,~] = FAIR(3,Prt_5i,f_45,'BE');
R.eta_tL = (ht_45-ht_5)/(ht_45-ht_5i);

%Mixer Calculations
R.A_prime = I.A/(lmbmelme2tlpf+I.e1+I.e2);
f_6 = f_45;
R.f_6A = f_6/(1+R.A_prime);
ht_16 = ht_13;
Prt_16 = Prt_13;
ht_6 = ht_5; 
ht_6A = (ht_6+R.A_prime *ht_16)/(1+R.A_prime);
R.tau_M = ht_6A/ht_6;
R.Pt16_Pt6 = I.pi_f/(I.pi_c*I.pi_b*R.pi_tH*R.pi_tL);
R.Tt_6 = R.Tt_5;

R.Tt_16 = R.Tt_13;
[Tt6_T6,Pt6_P6,R.MFP_6] = MASSFP_NR(R.Tt_6,f_45,I.M_6);
R.A45_A6 = R.pi_tL*sqrt(R.Tt_45/R.Tt_6)*(R.MFP_6/R.MFP_45);
PtP16 = Pt6_P6*R.Pt16_Pt6;
Pr_16 = Prt_16/PtP16;
[T_16,h_16,~,Phi_16,Cp_16,R_16,gamma_16,a_16] = FAIR(3,Pr_16,0,'BE');
V_16 = sqrt(2*g_c*(ht_16-h_16)*778.16);
R.M_16 = V_16/a_16;
[TtT16,PtP16,MFP_16] = MASSFP_NR(R.Tt_16,0,R.M_16);
R.A16_A6 = R.A_prime*sqrt(R.Tt_16/R.Tt_6)*(R.MFP_6/MFP_16)/R.Pt16_Pt6;
R.A6_A6A = 1/(1+R.A16_A6);
T_6 = R.Tt_6/Tt6_T6;
[~,h_6,Pr_6,Phi_6,Cp_6,R_6,gamma_6,a_6] = FAIR(1,T_6,f_6,'BE');
[R.M_6A,R.MFP_6A,R.Tt_6A,R.gamma_6A,R.Cp_6A] = MIXER_VSH(I.M_6,R.M_16,gamma_16,gamma_6,R_6,R.A_prime,R.A16_A6,T_6,R.f_6A,ht_6A);
R.pi_M_ideal = (1+R.A_prime)*sqrt(R.tau_M)*R.A6_A6A*R.MFP_6/R.MFP_6A;
R.pi_M = R.pi_M_ideal*I.pi_M_max;
R.f_AB = FUEL_NR(ht_6A,I.eta_AB,I.h_PR,I.Tt_7);

%Afterburner OFF
R.f_0_off = R.lmbme1me2*R.f/(1+I.A); % Overall fuel-air ratio (AB = off). [Dimensionless]
[~,ht_7_off,Prt_7_off,Phit_7_off,~,~,~,~] = FAIR(1,R.Tt_6A,R.f_0_off,'BE');
ht_9_off = ht_7_off;
Prt_9_off = Prt_7_off;
R.PtP_9_off = I.P0_P9*R.pi_r*R.pi_d*R.pi_cL*R.pi_cH*I.pi_b*R.pi_tH*R.pi_tL*R.pi_M*R.pi_ABdry*I.pi_n;
Pr_9_off = Prt_9_off/R.PtP_9_off;
[T_9_off,h_9_off,~,Phi_9_off,Cp_9_off,R_9_off,gamma_9_off,a_9_off] = FAIR(3,Pr_9_off,R.f_0_off,'BE');
V_9_off = sqrt(2*g_c*778.16*(ht_9_off-h_9_off));
R.M_9_off = V_9_off/a_9_off;
if R.M_9_off > 1
    R.M_8_off = 1;
else
    R.M_8_off = R.M_9_off;
end
[TtT_9_off,~,MFP_9_off] = MASSFP_NR(R.Tt_6A,R.f_6A,R.M_9_off);
[TtT_8_off,PtP_8_off,MFP_8_off] = MASSFP_NR(R.Tt_6A,R.f_6A,R.M_8_off);
R.A8_A6_dry = R.MFP_6*(1+R.A_prime)*sqrt(R.Tt_6A/R.Tt_6)/(R.pi_M*R.pi_ABdry*MFP_8_off);
R.A9_A8_off = MFP_8_off/(MFP_9_off*I.pi_n);
R.T9_T0_off =  T_9_off/R.T_0;
R.V9_V0_off = V_9_off/V_0;
R.M9_M0_off = R.M_9_off/I.M_0;
R.F_m_dot_off = (a_0/g_c)*((1+R.f_0_off - I.B/(1+I.A))*(V_9_off/a_0 + R_9_off*T_9_off*a_0*(1-I.P0_P9)/(R_0*R.T_0*V_9_off*gamma_0)) - I.M_0); % Uninstalled Specific Thrust (AB = on). [lbf/lbm*s]
R.S_off = 3600*R.f_0_off/R.F_m_dot_off; %Uninstalled TSFC (AB = on). [1/hr]
R.eta_P_off = (2*g_c*I.M_0*R.F_m_dot_off/a_0)/((1 + R.f_0_off - I.B/(1+I.A))*((V_9_off/a_0)^2)-I.M_0^2);
R.eta_TH_off =(((1+R.f_0_off-I.B/(1+I.A))*(V_9_off^2) - V_0^2)/(2*g_c*778.16)+(I.CTOL+I.CTOH)*R.h_0)/(R.f_0_off*I.h_PR);

%Afterburner ON
R.f_0_on = (R.lmbme1me2*R.f + (R.f_AB*(1+I.A-I.B)))/(1+I.A); % Overall fuel/air ratio (AB = on). [Dimensionless]
[~,ht_7,Prt_7,Phit_7,~,~,~,~] = FAIR(1,I.Tt_7,R.f_0_on,'BE');
ht_9 = ht_7;
Prt_9 = Prt_7;
R.PtP_9_on = I.P0_P9*R.pi_r*R.pi_d*R.pi_cL*R.pi_cH*I.pi_b*R.pi_tH*R.pi_tL*R.pi_M*I.pi_AB*I.pi_n;
Pr_9 = Prt_9/R.PtP_9_on;
[T_9,h_9,~,Phi_9,Cp_9,R_9,gamma_9,a_9] = FAIR(3,Pr_9,R.f_0_on,'BE');
V_9 = sqrt(2*g_c*778.16*(ht_9-h_9));
R.M_9_on = V_9/a_9;
if R.M_9_on > 1
    R.M_8_on = 1;
else
    R.M_8_on = R.M_9_on;
end
[TtT_9_on,~,MFP_9_on] = MASSFP_NR(I.Tt_7,R.f_0_on,R.M_9_on);
[TtT_8_on,PtP_8_on,MFP_8_on] = MASSFP_NR(I.Tt_7,R.f_0_on,R.M_8_on);
R.A9_A8_on = MFP_8_on/(MFP_9_on*I.pi_n);
R.T9_T0_on =  T_9/R.T_0;
R.V9_V0_on = V_9/V_0;
R.M9_M0_on = R.M_9_on/I.M_0;
R.F_m_dot_on = (a_0/g_c)*((1+R.f_0_on - I.B/(1+I.A))*(V_9/a_0 + R_9*T_9*a_0*(1-I.P0_P9)/(R_0*R.T_0*V_9*gamma_0)) - I.M_0); % Uninstalled Specific Thrust (AB = on). [lbf/lbm*s]
R.S_on = 3600*R.f_0_on/R.F_m_dot_on; %Uninstalled TSFC (AB = on). [1/hr]
R.eta_P_on = (2*g_c*I.M_0*R.F_m_dot_on/a_0)/((1 + R.f_0_on - I.B/(1+I.A))*((V_9/a_0)^2)-I.M_0^2);
R.eta_TH_on =(((1+R.f_0_on-I.B/(1+I.A))*(V_9^2) - V_0^2)/(2*g_c*778.16)+(I.CTOL+I.CTOH)*R.h_0)/(R.f_0_on*I.h_PR);
end
%% OFX Calculator
function O = OFX_LBTF_VSH(R,I,M_0,Alt,AB) %617
%Data Preallocations and Initial value setting.
O = struct('Alt',Alt,'A',0,'a_0',0,'f',R.f,'ht_0',0,'ht_2',0,'ht_3',0,'ht_4',0,'ht_45',0,...
    'M_0',M_0,'M_6',R.M_6,'M_6A',R.M_6A,'Cp_6A',0,'gamma_6A',0,...
    'm_dot',R.m_dot,'rho_0',0,...
    'P_0',0,'T_0',0,'Tt_4',0,'Tt_45',0,'Tt_5',0,'Pt16_Pt6',0,'pi_ABdry',0,'pi_c',R.pi_c*1.1,'pi_cH',R.pi_cH,...
    'pi_cL',R.pi_cL,'pi_d',0,'pi_f',R.pi_f,'pi_tH',R.pi_tH,...
    'pi_tL',R.pi_tL,'pi_r',0,'tau_cH',R.tau_cH,'tau_cL',R.tau_cL,'tau_f',R.tau_f,...
    'tau_m1',R.tau_m1,'tau_m2',R.tau_m2,'tau_r',0,'tau_tH',R.tau_tH,...
    'tau_tL',R.tau_tL,'theta_0',0,'V_0',0);

%Preliminary Computations
g_c = 32.174;
[O.T_0,O.a_0,O.P_0,O.rho_0] = ATMOSphere(O.Alt,'BE',I.day);
[~,h_0,Pr_0,Phi_0,Cp_0,R_0,gamma_0,a_0] = FAIR(1,O.T_0,0,'BE');
O.V_0 = O.M_0*a_0;
O.ht_0 = h_0 + (O.V_0^2)/(2*g_c*778.16);
O.tau_r = O.ht_0/h_0; %OFX Ram enthalpy ratio.
O.theta_0 = O.tau_r*O.T_0/518.67; %OFX Throttle Ratio.
[Tt_0,~,Prt_0,Phit_0,~,~,~,~] = FAIR(2,O.ht_0,0,'BE');
O.pi_r = Prt_0/Pr_0; %OFX Ram pressure ratio.
eta_R_spec = MIL_E_5008B(O.M_0);
O.pi_d = I.pi_d_max*eta_R_spec;
O.ht_2 = O.ht_0;
Prt_2 = Prt_0;
O.A = R.A*O.T_0*O.tau_r/(R.tau_r*R.T_0); %Initial guess for bypass.
O.pi_ABdry = 1-(1-I.pi_AB)/2;
if O.theta_0 >= R.theta_0
    Tt_4_guess = R.Tt_4;
else
    Tt_4_guess = R.Tt_4*O.theta_0/R.theta_0;
end
O.Tt_4 = Tt_4_guess;
[~,O.ht_4,Prt_4,Phit_4,~,~,~,~] = FAIR(1,O.Tt_4,O.f,'BE');
O.ht_45 = O.ht_4*O.tau_m1*O.tau_tH*O.tau_m2;
f_45 = O.f/(1+O.f+(R.e1+R.e2)/R.lmbme1me2);
[O.Tt_45,~,Prt_45,Phit_45,~,~,~,~] = FAIR(2,O.ht_45,f_45,'BE');
ht_5 = O.ht_45*O.tau_tL;
[O.Tt_5,~,Prt_5,Phit_5,~,~,~,~] = FAIR(2,ht_5,f_45,'BE');

%STAGE 1
k_max = 100;
k = 0;
k_throttle = 100;
while O.pi_c > R.pi_c && k_throttle > 0
    O.Tt_4 = Tt_4_guess*k_throttle/100;
    m_error = 1;
    while m_error > 0.001 && k <= k_max
        M6_error = 1;
        while M6_error > 0.0005 && k <= k_max
            a_error = 1;
            while a_error > 0.001 && k <= k_max
                req3 = true;
                while req3 == true && k <= k_max
                    O.ht_3 = O.ht_0*O.tau_cL*O.tau_cH;
                    [Tt_3,~,Prt_3,Phit_3,~,~,~,~] = FAIR(2,O.ht_3,0,'BE');
                    O.A_prime = O.A/((1+O.f)*R.lmbme1me2+R.e1+R.e2);%%%%%%%%%%
                    [~,O.ht_4,Prt_4,Phit_4,~,~,~,~] = FAIR(1,O.Tt_4,O.f,'BE');
                    [O.pi_tH,O.tau_tH,O.Tt_45] = TURBC(O.Tt_4,O.ht_4,O.f,R.A4_A45,R.M_4,R.M_45,R.eta_tH,O.Tt_45,O.ht_3,R.lmbme1me2,R.e1,R.e2);
                    [O.pi_tL,O.tau_tL,O.Tt_5] = TURB(O.Tt_45,f_45,R.A45_A6,R.M_45,O.M_6,R.eta_tL,O.Tt_5);
                    [~,ht_5,Prt_5,Phit_5,~,~,~,~] = FAIR(1,O.Tt_5,f_45,'BE');
                    O.tau_lambda = O.ht_4/h_0;
                    O.tau_f = 1 + ((1-O.tau_tL)*I.eta_mL*(R.lmbme1me2*(1+O.f)*O.tau_lambda*O.tau_tH/O.tau_r+(R.e1*O.tau_tH+R.e2)*O.tau_cL*O.tau_cH)-(1+O.A)*R.P_TOL/(O.tau_r*I.eta_mPL*O.m_dot*h_0))/((R.tau_cL-1)/(R.tau_f-1)+O.A);
                    O.tau_cL = 1 + (O.tau_f-1)*(R.tau_cL-1)/(R.tau_f-1);
                    O.tau_cH = (1 + (1-O.tau_tH)*I.eta_mH*(R.lmbme1me2*(1+O.f)*O.tau_lambda/(O.tau_r*O.tau_cL)) - (1+O.A)*R.P_TOH/(O.tau_r*O.tau_cL*I.eta_mPH*O.m_dot*h_0))/(1-R.e1*(1-O.tau_tH)*I.eta_mH);
                    ht_13 = O.ht_2*O.tau_f;
                    ht_25 = O.ht_2*O.tau_cL;
                    O.ht_3 = ht_25*O.tau_cH;
                    ht_13i = O.ht_2*(1+R.eta_f*(O.tau_f-1));
                    ht_25i = O.ht_2*(1+R.eta_cL*(O.tau_cL-1));
                    ht_3i = ht_25*(1+R.eta_cH*(O.tau_cH-1));
                    [Tt_13,~,Prt_13,Phit_13,~,~,~,~] = FAIR(2,ht_13,0,'BE');
                    [Tt_25,~,Prt_25,Phit_25,~,~,~,~] = FAIR(2,ht_25,0,'BE');
                    [Tt_3,~,Prt_3,Phit_3,~,~,~,~] = FAIR(2,O.ht_3,0,'BE');
                    [Tt_13i,~,Prt_13i,Phit_13i,~,~,~,~] = FAIR(2,ht_13i,0,'BE');
                    [Tt_25i,~,Prt_25i,Phit_25i,~,~,~,~] = FAIR(2,ht_25i,0,'BE');
                    [Tt_3i,~,Prt_3i,Phit_3i,~,~,~,~] = FAIR(2,ht_3i,0,'BE');
                    O.pi_f = Prt_13i/Prt_2;
                    O.pi_cL = Prt_25i/Prt_2;
                    O.pi_cH = Prt_3i/Prt_25;
                    O.pi_c = O.pi_cL*O.pi_cH;
                    O.f = FUEL_NR(O.ht_3,I.eta_b,I.h_PR,O.Tt_4);
                    O.tau_m1 = (R.lmbme1me2*(1+O.f)+R.e1*O.tau_r*O.tau_cL*O.tau_cH/O.tau_lambda)/(R.lmbme1me2*(1+O.f)+R.e1);
                    O.tau_m2 = (R.lmbme1me2*(1+O.f)+R.e1+R.e2*(O.tau_r*O.tau_cL*O.tau_cH/(O.tau_lambda*O.tau_m1*O.tau_tH)))/(R.lmbme1me2*(1+O.f)+R.e1+R.e2);
                    ht_6 = ht_5;
                    Tt_6 = O.Tt_5;
                    ht_16 = ht_13;
                    Tt_16 = Tt_13;
                    Pt_6 = O.P_0*O.pi_r*O.pi_d*O.pi_cL*O.pi_cH*I.pi_b*O.pi_tH*O.pi_tL;
                    f_45 = O.f/(1+O.f+(R.e1+R.e2)/R.lmbme1me2);
                    [TtT6,PtP6,O.MFP_6] = MASSFP_NR(Tt_6,f_45,O.M_6);
                    P_6 = Pt_6/PtP6;
                    T_6 = Tt_6/TtT6;
                    Pt_16 = O.P_0*O.pi_r*O.pi_d*O.pi_f;
                    PtP16 = Pt_16/P_6; %Kutta condition
                    [~,PtP,~] = MASSFP_NR(Tt_16,0,1); %PtP at M = 1
                    if PtP16 > PtP || PtP16 < 1
                        O.M_6 = O.M_6-0.01;
                        k = k+1;
                    else
                        req3 = false;
                    end
                    %fprintf('loop1\n')
                end
                [O.M_16,TtT16,~,MFP_16] = RGCOMPR_NR(3,Tt_16,0,PtP16);
                T_16 = Tt_16/TtT16;
                A_prime = R.A16_A6*Pt_16*MFP_16*sqrt(Tt_6/Tt_16)/(Pt_6*O.MFP_6);
                O.Pt16_Pt6 = Pt_16/Pt_6;
                a_error = abs((A_prime-O.A_prime)/O.A_prime);
                O.A = A_prime*((1+O.f)*R.lmbme1me2+R.e1+R.e2);
                k = k+1;
                %fprintf('loop2\n')
            end
            [~,h_6,Pr_6,Phi_6,Cp_6,R_6,gamma_6,a_6] = FAIR(1,T_6,f_45,'BE');
            [~,h_16,Pr_16,Phi_16,Cp_16,R_16,gamma_16,a_16] = FAIR(1,T_16,0,'BE');
            ht_6A = (ht_6+O.A_prime*ht_16)/(1+O.A_prime);
            O.tau_M = ht_6A/ht_6;
            O.f_6A = f_45/(1+O.A_prime);
            [O.M_6A,MFP_6A,Tt_6A,O.gamma_6A,O.Cp_6A] = MIXER_VSH(O.M_6,O.M_16,gamma_16,gamma_6,R_6,O.A_prime,R.A16_A6,T_6,O.f_6A,ht_6A);
            O.pi_M_ideal = sqrt(Tt_6A/Tt_6)*O.MFP_6*(1+O.A_prime)/(MFP_6A*(1+R.A16_A6));
            O.pi_M = O.pi_M_ideal*I.pi_M_max;
            PtP9 = O.pi_r*O.pi_d*O.pi_cL*O.pi_cH*I.pi_b*O.pi_tH*O.pi_tL*O.pi_M*O.pi_ABdry*I.pi_n*I.P0_P9;
            [O.M_9,TtT9,~,MFP_9] = RGCOMPR_NR(3,Tt_6A,O.f_6A,PtP9);
            if O.M_9 > 1
                O.M_8 = 1;
            else
                O.M_8 = O.M_9;
            end
            [TtT_8,PtP_8,MFP_8] = MASSFP_NR(Tt_6A,O.f_6A,O.M_8);
            O.MFP_6 =  MFP_8*O.pi_M*O.pi_ABdry*R.A8_A6_dry*sqrt(Tt_6/Tt_6A)/(1+O.A_prime);
            [M_6new,TtT_6,PtP_6,~] = RGCOMPR_NR(5,Tt_6,f_45,O.MFP_6);
            M6_error = abs(M_6new-O.M_6);
            if M6_error > 0.0005
                if O.M_6 > M_6new
                    O.M_6 = O.M_6 - 0.0001;
                else
                    O.M_6 = O.M_6 + 0.002;
                end
            end
            %fprintf('loop3\n')
        end
        [~,TtT,PtP,MFP_4] = RGCOMPR_NR(1,O.Tt_4,O.f,1);
        m_new = R.m_dot*(1+R.f)*O.P_0*(1+O.A)*O.pi_r*O.pi_d*O.pi_c*MFP_4*sqrt(R.Tt_4/O.Tt_4)/((1+O.f)*(R.P_0*(1+R.A)*R.pi_r*R.pi_d*R.pi_c*R.MFP_4));
        m_error = abs((m_new-O.m_dot)/O.m_dot);
        O.m_dot = m_new;
        k = k + 1;
        %fprintf('loop4\n')
    end
    k_throttle = k_throttle - 1;
end
    O.Tt_7 = AB*(R.Tt_7-Tt_6A) + Tt_6A;
    if AB == 0
        O.f_AB =0;
    else
        O.f_AB = FUEL_NR(ht_6A,I.eta_b,I.h_PR,O.Tt_7);
    end
    O.f_7 = O.f_6A + O.f_AB;
    O.pi_AB = O.pi_ABdry + 0.01*AB*(I.pi_AB-O.pi_ABdry);
    PtP9 = O.pi_r*O.pi_d*O.pi_cL*O.pi_cH*I.pi_b*O.pi_tH*O.pi_tL*O.pi_M*O.pi_AB*I.pi_n*I.P0_P9;
    Tt_9 = O.Tt_7;
    [M_9,TtT9,~,MFP_9] = RGCOMPR_NR(3,Tt_9,O.f_7,PtP9);
    m_dot_9 = O.m_dot*(1+O.f_7)*(1-R.B/(1+O.A));
    Pt_9 = O.P_0*O.pi_r*O.pi_d*O.pi_cL*O.pi_cH*I.pi_b*O.pi_tH*O.pi_tL*O.pi_M*O.pi_AB*I.pi_n;
    O.A_9 = m_dot_9*sqrt(Tt_9)/(Pt_9*MFP_9);
    T_9 = Tt_9/TtT9;
    [~,h_9,Pr_9,Phi_9,Cp_9,R_9,gamma_9,a_9] = FAIR(1,T_9,O.f_7,'BE');
    V_9 = M_9*a_9;
    O.f_0 = (O.f*R.lmbme1me2 + O.f_AB*(1+O.A-R.B))/(1+O.A);
    O.F_m_dot = (a_0/g_c)*((1+O.f_0 - I.B/(1+O.A))*(V_9/a_0 + R_9*T_9*a_0*(1-I.P0_P9)/(R_0*R.T_0*V_9*gamma_0)) - O.M_0); % Uninstalled Specific Thrust (AB = on). [lbf/lbm*s]
    O.S = 3600*O.f_0/O.F_m_dot;
    O.F = O.m_dot*O.F_m_dot;
    [~,TtT0,PtP0,MFP_0] = RGCOMPR_NR(1,Tt_0,0,O.M_0);
    O.A_0 = O.m_dot*sqrt(Tt_0)/(O.P_0*PtP0*MFP_0);
    O.LP_RPM = 100*sqrt(h_0*O.tau_r*(O.tau_f-1)/(R.h_0*R.tau_r*(R.tau_f-1)));
    O.HP_RPM = 100*sqrt(h_0*O.tau_r*O.tau_cL*(O.tau_cH-1)/(R.h_0*R.tau_r*R.tau_cL*(R.tau_cH-1)));
    O.eta_P = (2*g_c*O.M_0*O.F_m_dot/a_0)/((1 + O.f_0 - R.B/(1+O.A))*((V_9/a_0)^2)-O.M_0^2);
    O.eta_TH =(((1+O.f_0-I.B/(1+O.A))*(V_9^2) - O.V_0^2)/(2*g_c*778.16)+(I.CTOL+I.CTOH)*h_0)/(O.f_0*I.h_PR);
    %fprintf('success!\n')
    %fprintf('%5f',k_throttle)
end
%% ATMOS 
%TO-DO:
%1) Code-in the cold day model.
%2) Code-in the tropical day model.

%----------------------------------------
%AE 435: Air-Breathing Preliminary Design
%RFP: AEDsys "ATMOS" Code
%By Jesus Ferrand
%Submitted to: Mark Ricklick Phd.
%----------------------------------------

%NOTES:
%AEDsys' ATMOS program only has data up to h = 30km (98.4252 kft) for its
%standard day model. The data for the cold, hot, and tropic days is limited
%to heights of 22.5km (73.8189 kft), 20.5km (67.2572 kft), and 21km
%(68.8976 kft) respectively. This code will accept inputs outside these
%height limits but bear in mind that the output will be extrapolated then.

%CODE:
%==========================================================================
%[T,a,P,rho] = ATMOSphere(10000       ,'SI','Standard');
%[T,a,P,rho] = ATMOSphere(10000/0.3048,'BE','Standard');

function [T,a,P,rho] = ATMOSphere(h,units,day)
if strcmpi(units,'BE') == 1
    h = h*0.3048; % Input height in BE units converted to SI units. [m]
end
r_0 = 6356.577*1000; % Earth's radius. [m]
g_0 = 9.80665; % Earth's acceleration due to gravity at its surface. [m/s^2]
R = 8.31432; % Universal gas constant. [J/mol-K]
W_0 = 28.9644/1000; % Air molecular weight. [kg/mol]
P_0 = 101325; % SLS pressure. [Pa]
z = h*r_0/(r_0+h); % "Height Parameter." [m]

if strcmpi(day,'Standard') == 1
    T_0 = 288.15; % Surface Temperature. [K]
    c = [0*1000,-6.5/1000;...
        11*1000,0/1000;...
        20*1000,1/1000;...
        32*1000,2.8/1000;...
        47*1000,0/1000;...
        51*1000,-2.8/1000;...
        71*1000,-2/1000;...
        84.582*1000,0/1000];
    
    T_1 = T_0 + c(1,2)*(c(2,1)-c(1,1)); % [K]
    T_2 = T_1 + c(2,2)*(c(3,1)-c(2,1)); % [K]
    T_3 = T_2 + c(3,2)*(c(4,1)-c(3,1)); % [K]
    T_4 = T_3 + c(4,2)*(c(5,1)-c(4,1)); % [K]
    T_5 = T_4 + c(5,2)*(c(6,1)-c(5,1)); % [K]
    T_6 = T_5 + c(6,2)*(c(7,1)-c(6,1)); % [K]
    
    P_1 = P_0*power(T_0/T_1,g_0*W_0/(R*c(1,2))); % [Pa]
    P_2 = P_1*exp(-g_0*W_0*(c(3,1)-c(2,1))/(R*T_1)); % [Pa]
    P_3 = P_2*power(T_2/T_3,g_0*W_0/(R*c(3,2))); % [Pa]
    P_4 = P_3*power(T_3/T_4,g_0*W_0/(R*c(4,2))); % [Pa]
    P_5 = P_4*exp(-g_0*W_0*(c(6,1)-c(5,1))/(R*T_4)); % [Pa]
    P_6 = P_5*power(T_5/T_6,g_0*W_0/(R*c(6,2))); % [Pa]
    %P_7 = P_6*power(T_6/T_7,g_0*W_0/(R*c(7,2))); % [Pa]
    
    if z >= 0 && z < c(2,1)
        T_i = T_0;
        P_i = P_0;
        L_i = c(1,2);
        z_i = c(1,1);
    elseif z >= c(2,1) && z < c(3,1)
        T_i = T_1;
        P_i = P_1;
        L_i = c(2,2);
        z_i = c(2,1);
    elseif z >= c(3,1) && z < c(4,1)
        T_i = T_2;
        P_i = P_2;
        L_i = c(3,2);
        z_i = c(3,1);
    elseif z >= c(4,1) && z < c(5,1)
        T_i = T_3;
        P_i = P_3;
        L_i = c(4,2);
        z_i = c(4,1);
    elseif z >= c(5,1) && z < c(6,1)
        T_i = T_4;
        P_i = P_4;
        L_i = c(5,2);
        z_i = c(5,1);
    elseif z >= c(6,1) && z < c(7,1)
        T_i = T_5;
        P_i = P_5;
        L_i = c(6,2);
        z_i = c(6,1);
    elseif z >= c(7,1) && z < c(8,1)
        T_i = T_6;
        P_i = P_6;
        L_i = c(7,2);
        z_i = c(7,1);
    end
    T = T_i + L_i*(z-z_i); % Temperature @ input height. [K]
    P = pressure_lapse(L_i,T,T_i,P_i,z,z_i,R,g_0,W_0); % Static Pressure @ input height. [Pa]
    rho = P*W_0/(R*T); % Static Density @ input height. [kg/m^3]
    a = sqrt(1.4*R*T/W_0); % Speed of Sound @ input height. [m/s]
    
elseif strcmpi(day,'Hot') == 1
    T_0 = 312.60; % Surface Temperature. [K]
    c = [0*1000,-7/1000;...
        12*1000,0.8/1000;...
        20.5*1000,1.4/1000];
    T_1 = T_0 + c(1,2)*(c(2,1)-c(1,1)); % [K]
    T_2 = T_1 + c(2,2)*(c(3,1)-c(2,1)); % [K]
    P_1 = P_0*power(T_0/T_1,g_0*W_0/(R*c(1,2))); % [Pa]
    P_2 = P_1*power(T_1/T_2,g_0*W_0/(R*c(2,2))); % [Pa]
    if h >= 0 && h < c(2,1) %0km <= h < 12km
        T_i = T_0;
        P_i = P_0;
        L_i = c(1,2);
        h_i = c(1,1);
    elseif h >= c(2,1) && h < c(3,1) %12km <= h < 20.5km
        T_i = T_1;
        P_i = P_1;
        L_i = c(2,2);
        h_i = c(2,1);
    elseif h >= c(3,1) %20.5km <= h
        T_i = T_2;
        P_i = P_2;
        L_i = c(3,2);
        h_i = c(3,1);
    end
    z_i = h_i*r_0/(r_0+h_i);
    T = T_i + L_i*(h-h_i); % Temperature @ input height. [K]
    P = pressure_lapse(L_i,T,T_i,P_i,z,z_i,R,g_0,W_0); % Static Pressure @ input height. [Pa]
    rho = P*W_0/(R*T); % Static Density @ input height. [kg/m^3]
    a = sqrt(1.4*R*T/W_0); % Speed of Sound @ input height. [m/s]
elseif strcmpi(day,'Cold') == 1
    T_0 = 222.1; % Surface Temperature. [K]
elseif strcmpi(day,'Tropical') == 1
    T_0 = 305.27; % Surface Temperature. [K]
end
if strcmpi(units,'BE') == 1
    T = T*1.8; % Static Temperature @ input height. [R]
    a = a/0.3048; % Speed of sound @ input height. [ft/s]
    P = P*2116.2/101325; % Static Pressure @ input height. [psf]
    rho = rho/((100^3)/((12*2.54)^3))/(0.45359*32.174); % Static Density @ input height. [slug/ft^3]
end
    function P = pressure_lapse(L_i,T,T_i,P_i,z,z_i,R,g_0,W_0)
        if L_i == 0
            P = P_i*exp(-g_0*W_0*(z-z_i)/(R*T_i));
        else
            P = P_i*power(T_i/T,g_0*W_0/(R*L_i));
        end
    end
end
%% Drag Model
function [K1_future,K1_current,CD0_current,CD0_future] = K1_CD0(M)

%calculate K1 based on graph on page.37 of text
if M<=1
    K1_future = 0.18;
elseif M>1
    K1_future = 0.18*M;
end

if M<=0.8
    K1_current = 0.14;
elseif M<=1.2
    K1_current = 0.15*M + 0.02;
elseif M>1.2
    K1_current = 0.375*M - 0.25;
end

%calculate CD0 max and mil from graph on page 37
if M<=0.8
    CD0_current = 0.018;
elseif M<=1.2
    CD0_current = 0.055*M - 0.026;
elseif M>1.2
    CD0_current = -0.0025*M + 0.043;
end

if M<=0.8
    CD0_future = 0.014;
elseif M<=0.9
    CD0_future = 0.02*M - 0.002;
elseif M<=1.1
    CD0_future = 0.05*M - 0.029;
elseif M<=1.2
    CD0_future = 0.02*M + 0.004;
elseif M>1.2
    CD0_future = 0.028;
end
end
%% FAIR: Variable Specific Heat (VSH) Model 
%----------------------------------------------
%Title: "Mattingly's FAIR subroutine in MATLAB"
%By: Jesus Ferrand
%Credits: Jack D. Mattingly
%----------------------------------------------
%INFO: MATLAB version of Mattingly's FAIR subroutine as explained in his
%work "Elements of Gas Turbine Propulsion." This subroutine (now MATLAB
%function) resolves the thermodynamic state of a gas mixture of air and
%combustion products between air and fuels of the family "(CH)n" with
%variable gas properties. It takes one of four thermodynamic inputs:
%Temperature(T), Enthalpy (h), specific thermal entropy (Phi), or Reduced
%Pressure (Pr).

%INPUTS:
%   mode: INTEGER values 1,2,3, and 4!
%       mode 1: Temperature Known
%       mode 2: Enthalpy Known
%       mode 3: Reduced Pressure Known
%       mode 4: Specific Entropy Known
%   IN: POSITIVE, FLOAT, values denoting the thermodynamic input.
%   f: POSITIVE, FLOAT, values denoting the fuel-to-air ratio.
%   units: STRING either 'BE' or 'SI' to specify units.
%OUTPUTS:
%   [T,h,Pr,Phi,Cp,R,gamma,a]: ARRAY
%   T: Temperature. [R or K]
%   h: Enthalpy. [Btu/lbm or J/kg]
%   Pr: Reduced Pressure. [Dimensionless]
%   Phi: Specific Entropy. [Btu/lbmR or J/kgK]
%   Cp: Specific Heat Capacity at constant pressure. [Btu/lbm or J/kg]
%   R: Specific Gas constant. [Btu/lbm or J/kg]
%   gamma: Ratio of Specific Heats. (gamma = Cp/Cv) [Dimensionless]
%   a: Local Speed of Sound. [ft/s or m/s]

function [T,h,Pr,Phi,Cp,R,gamma,a] = FAIR(mode,IN,f,units)
if strcmpi(units,'SI')
    T_factor = 1.8; %[R]/[K]
    h_factor = 1/(2326); %[Btu/lbm]/[J/kg]
    Cp_factor = 1/(4186.8); %[Btu/lbm-R]/[J/kg-K]
    a_factor = 0.3048; %[m]/[ft]
else
    T_factor = 1;
    h_factor = 1;
    Cp_factor = 1;
    a_factor = 1;
end
C0 = (f*+7.3816638e-2 +2.5020051e-1)/(1+f);
C1 = (f*+1.2258630e-3 -5.1536879e-5)/(1+f); D1 = C1/2;
C2 = (f*-1.3771901e-6 +6.5519486e-8)/(1+f); D2 = C2/3; S2 = C2/2;
C3 = (f*+9.9686793e-10 -6.7178376e-12)/(1+f); D3 = C3/4; S3 = C3/3;
C4 = (f*-4.2051104e-13 -1.5128259e-14)/(1+f); D4 = C4/5; S4 = C4/4;
C5 = (f*+1.0212913e-16 +7.6215767e-18)/(1+f); D5 = C5/6; S5 = C5/5;
C6 = (f*-1.3335668e-20 -1.4526770e-21)/(1+f); D6 = C6/7; S6 = C6/6;
C7 = (f*+7.2678710e-25 +1.0115540e-25)/(1+f); D7 = C7/8; S7 = C7/7;
h_R= (f*30.58153 -1.7558886)/(1+f);
Phi_ref = (f*0.6483398+0.0454323)/(1+f);
R = 1.9857117/(28.97-f*0.946186);
Phi_0 = 1.5784416522122042; %Reference Phi used in Reduced Pressure.
diff = 0.5;
k = 0; % Counter variable for numeric solver convergence criterion.
k_max = 100; % Maximum number of iterations.
switch mode
    case 1 % Temperature Known
        T = IN*T_factor ;
        R= R/Cp_factor;
        Cp = (C0 + C1*T + C2*T^2 + C3*T^3 + C4*T^4 + C5*T^5 + C6*T^6 + C7*T^7)/Cp_factor;
        h = (h_R + C0*T + D1*T^2 + D2*T^3 + D3*T^4 + D4*T^5 + D5*T^6 + D6*T^7 + D7*T^8)/h_factor;
        Phi = (Phi_ref + log(T)*C0 + C1*T + S2*T^2 + S3*T^3 + S4*T^4 + S5*T^5 + S6*T^6 + S7*T^7)/Cp_factor;
        Pr = exp((Phi-Phi_0/Cp_factor)/R);
        gamma = Cp/(Cp-R);
        a = sqrt(gamma*R*T*32.174*778.16)*a_factor;
        T = T/T_factor;
    case 2 % Enthalpy Known
        T = IN/0.26;
        h = IN*h_factor;
        h_0 = h_R - h;
        while diff > 0.00000000000001 && k <= k_max
            T0 = T;
            h_iter = h_0 + T0*C0 + D1*T0^2 + D2*T0^3 + D3*T0^4 + D4*T0^5 + D5*T0^6 + D6*T0^7 + D7*T0^8;
            Cp = C0 + T0*C1 + C2*T0^2 + C3*T0^3 + C4*T0^4 + C5*T0^5 + C6*T0^6 + C7*T0^7;
            T = T0 - h_iter/Cp;
            diff = abs(T - T0);
            k = k + 1;
        end
        gamma = Cp/(Cp-R);
        Phi = (Phi_ref + log(T)*C0 + C1*T + S2*T^2 + S3*T^3 + S4*T^4 + S5*T^5 + S6*T^6 + S7*T^7)/Cp_factor;
        Pr = exp((Phi-Phi_0/Cp_factor)/R);
        a = sqrt(gamma*R*T*32.174*778.16)*a_factor;
    case {3,4} % Reduced Pressure Known (3) or Entropy Known (4)
        T = 1200; % Default Initial Guess for modes 2,3, and 4. [R]
        if mode == 3
            Pr = IN;
            Phi_C = Phi_ref-(log(IN)*R + Phi_0);
        else %mode == 4
            Phi_C = Phi_ref-(IN*Cp_factor);
        end
        while diff > 0.00000000001 && k <= k_max
            T0 = T;
            Phi = Phi_C + log(T0)*C0 + C1*T0 + S2*T0^2 + S3*T0^3 + S4*T0^4 + S5*T0^5 + S6*T0^6 + S7*T0^7;
            Phi_der = C0/T0 + C1 + C2*T0 + C3*T0^2 + C4*T0^3 + C5*T0^4 + C6*T0^5 + C7*T0^6;
            T = T0 - Phi/Phi_der;
            diff = abs(T - T0);
            k = k + 1;
        end
        R= R/Cp_factor;
        Cp = (C0 + C1*T + C2*T^2 + C3*T^3 + C4*T^4 + C5*T^5 + C6*T^6 + C7*T^7)/Cp_factor;
        h = (h_R + C0*T + D1*T^2 + D2*T^3 + D3*T^4 + D4*T^5 + D5*T^6 + D6*T^7 + D7*T^8)/h_factor;
        gamma = Cp/(Cp-R);
        a = sqrt(gamma*R*T*32.174*778.16)*a_factor;
        T = T/T_factor;
end
end
%% FUEL_NR: Fuel-Air-Ratio Numeric Solver
function f = FUEL_NR(enthalpy,eta_b,h_PR,Tt)
A0 = +2.5020051e-1;B0 = +7.3816638e-2;
A1 = -5.1536879e-5;B1 = +1.2258630e-3;
A2 = +6.5519486e-8;B2 = -1.3771901e-6;
A3 = -6.7178376e-12;B3 = +9.9686793e-10;
A4 = -1.5128259e-14;B4 = -4.2051104e-13;
A5 = +7.6215767e-18;B5 = +1.0212913e-16;
A6 = -1.4526770e-21;B6 = -1.3335668e-20;
A7 = +1.0115540e-25;B7 = +7.2678710e-25;
Ah_ref = -1.7558886;Bh_ref = 30.58153; % [Btu/lbm]

k = 0;
k_max = 100;
error = 0.5;

Asum = sum([Ah_ref,A0*Tt,A1*(Tt^2)/2,A2*(Tt^3)/3,A3*(Tt^4)/4,A4*(Tt^5)/5,A5*(Tt^6)/6,A6*(Tt^7)/7,A7*(Tt^8)/8]);
Bsum = sum([Bh_ref,B0*Tt,B1*(Tt^2)/2,B2*(Tt^3)/3,B3*(Tt^4)/4,B4*(Tt^5)/5,B5*(Tt^6)/6,B6*(Tt^7)/7,B7*(Tt^8)/8]);
AMH = Asum - enthalpy; %"Asum MINUS enthalpy"
BMH = Bsum - enthalpy; %"Bsum MINUS enthalpy"
FMB = eta_b*h_PR - Bsum; %"Fuel - Bsum"
FMA = eta_b*h_PR - Asum; %"Fuel - Asum"
f = (Tt*0.3-enthalpy)/(eta_b*h_PR - Tt*0.3); %Initial Guess
while error > 0.0000000000000001 && k <= k_max
    num = AMH + f*BMH;
    den = FMA + f*FMB;
    f_iter = num/den - f;
    f_prime = (BMH*den-num*FMB)/(den^2) - 1;
    f_new = f - f_iter/f_prime;
    error = abs(f_new-f);
    k = k+1;
    f = f_new;
end
end
%% MASSFP_NR: Mass Flow Parameter (MFP) Numeric Solver
function [TtT,PtP,MFP] = MASSFP_NR(Tt,f,M)
C0 = (f*+7.3816638e-2 +2.5020051e-1)/(1+f);
C1 = (f*+1.2258630e-3 -5.1536879e-5)/(1+f); D1 = C1/2;
C2 = (f*-1.3771901e-6 +6.5519486e-8)/(1+f); D2 = C2/3;
C3 = (f*+9.9686793e-10 -6.7178376e-12)/(1+f); D3 = C3/4; T3 = C3*2;
C4 = (f*-4.2051104e-13 -1.5128259e-14)/(1+f); D4 = C4/5; T4 = C4*3;
C5 = (f*+1.0212913e-16 +7.6215767e-18)/(1+f); D5 = C5/6; T5 = C5*4;
C6 = (f*-1.3335668e-20 -1.4526770e-21)/(1+f); D6 = C6/7; T6 = C6*5;
C7 = (f*+7.2678710e-25 +1.0115540e-25)/(1+f); D7 = C7/8; T7 = C7*6;
h_ref = (f*30.58153 -1.7558886)/(1+f); 
[~,ht,Prt,~,~,R,~,~] = FAIR(1,Tt,f,'BE');
K1 = (M^2)*R/2;
k = 0;
k_max = 100;
diff = 0.5;
if M > 1
    T = Tt*0.7;
else
    T = Tt*0.95;
end
while diff > 0.000001 && k < k_max
    Cp = C0 + C1*T + C2*(T^2) + C3*(T^3) + C4*(T^4) + C5*(T^5) + C6*(T^6)...
        + C7*(T^7);
    h_iter = h_ref + C0*T + D1*(T^2) + D2*(T^3) + D3*(T^4) + D4*(T^5)...
        + D5*(T^6) + D6*(T^7) + D7*(T^8) - ht + K1*Cp*T/(Cp-R);
    h_prime = Cp + K1*(1+(R*(C0 - C2*(T^2) - T3*(T^3) - T4*(T^4)...
        - T5*(T^5) - T6*(T^6) - T7*(T^7))-R^2)/power(Cp-R,2));
    Tnew = T - h_iter/h_prime;
    k = k + 1;
    diff = abs(Tnew-T);
    T = Tnew;
end
TtT = Tt/T;
[~,~,Pr,~,~,~,gamma,~] = FAIR(1,T,f,'BE');
PtP = Prt/Pr;
MFP = (M/PtP)*sqrt(gamma*TtT*32.174/R);
end
%% INV_MASSP_NR: Inverse MFP Solver
function [TtT,PtP,M] = INV_MASSP_NR(mode,Tt,f,MFP)
if mode == 4 %Supersonic solution Desired
    T = Tt*0.6;
else %Subsonic solution Desired.
    T = Tt*0.95;
end
k = 0;
k_max = 100;
diff = 0.5;
%Pre-computed values:
[~,ht,Prt,~,~,R,~,~] = FAIR(1,Tt,f,'BE');
C0 = (f*+7.3816638e-2 +2.5020051e-1)/(1+f); 
C1 = (f*+1.2258630e-3 -5.1536879e-5)/(1+f); D1 = C1/2;
C2 = (f*-1.3771901e-6 +6.5519486e-8)/(1+f); D2 = C2/3; S2 = C2/2; 
C3 = (f*+9.9686793e-10 -6.7178376e-12)/(1+f); D3 = C3/4; S3 = C3/3;  
C4 = (f*-4.2051104e-13 -1.5128259e-14)/(1+f); D4 = C4/5; S4 = C4/4; 
C5 = (f*+1.0212913e-16 +7.6215767e-18)/(1+f); D5 = C5/6; S5 = C5/5; 
C6 = (f*-1.3335668e-20 -1.4526770e-21)/(1+f); D6 = C6/7; S6 = C6/6; 
C7 = (f*+7.2678710e-25 +1.0115540e-25)/(1+f); D7 = C7/8; S7 = C7/7; 
h_R = (f*30.58153 -1.7558886)/(1+f);
Phi_R = (f*0.6483398+0.0454323)/(1+f);
Phi_R0 =  Phi_R - 1.5784416522122042;
C = R*(MFP*Prt)/sqrt(64.348*Tt);
while diff > 0.000001 && k < k_max    
    Cp = C7*T^7 + C6*T^6 + C5*T^5 + C4*T^4 + C3*T^3 + C2*T^2 + C1*T + C0;
    Pr = exp((Phi_R0 + C0*log(T) + S2*T^2 + S3*T^3 + S4*T^4 + S5*T^5 + S6*T^6 + S7*T^7 + C1*T)/R);
    htmh = -(D7*T^8 + D6*T^7 + D5*T^6 + D4*T^5 + D3*T^4 + D2*T^3 + D1*T^2 + C0*T + h_R - ht);
    dhdT2 = (htmh/T^2)^(1/2);
    f = Pr*dhdT2 - C;
    f_der = Pr*((dhdT2*(C0/T + C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4 + C6*T^5 + C7*T^6))/R - ((Cp/T^2 + (2*htmh)/T^3))/(2*dhdT2));
    Tn = T-f/f_der;
    diff = abs(Tn-T);
    T = Tn;
    k = k+1;
end
TtT = Tt/T;
PtP = Prt/Pr;
M = sqrt(2*(htmh)/(R*T*Cp/(Cp-R)));
end
%% TURB: Uncooled Turbine Solver
function [pi_tL,tau_tL,Tt_e] = TURB(Tt_i,f,Ai_Ae,M_i,M_e,eta_tL,Tt_eR)
[~,ht_i,Prt_i,Phit_i,~,~,~,~] = FAIR(1,Tt_i,f,'BE');
[TtTi,PtPi,MFP_i] = MASSFP_NR(Tt_i,f,M_i);
Tt_e = Tt_eR;
k = 0;
k_max = 100;
diff = 0.5;
while diff > 0.01 && k < k_max
    [TtTe,PtPe,MFP_e] = MASSFP_NR(Tt_e,f,M_e);
    pi_tL = MFP_i*Ai_Ae*sqrt(Tt_e/Tt_i)/MFP_e;
    Prt_ei = pi_tL*Prt_i;
    [~,ht_ei,Pr_ei,Phi_ei,~,~,~,~] = FAIR(3,Prt_ei,f,'BE');
    ht_e = ht_i - eta_tL*(ht_i-ht_ei);
    tau_tL = ht_e/ht_i;
    [Tt_en,~,~,~,~,~,~,~] = FAIR(2,ht_e,f,'BE');
    diff = abs(Tt_en-Tt_e);
    Tt_e = Tt_en;
    k = k+1;
end
end
%% TURBC: Cooled Turbine Solver
function [pi_tH,tau_tH,Tt_45] = TURBC(Tt_4,ht_4,f,A4_A45,M_4,M_45,eta_tH,Tt_45R,ht_3,lmbme1me2,e1,e2)
[TtT4,PtP4,MFP_4] = MASSFP_NR(Tt_4,f,M_4);
m_f = f*lmbme1me2;
m_4 = (1+f)*lmbme1me2;
m_41 = m_4+e1;
m_45 = m_41+e2;
f_41 = m_f/(m_41-m_f);
f_45 = m_f/(m_45-m_f);
ht_41 = (m_4*ht_4+e1*ht_3)/(m_4+e1);
[Tt_41,~,Prt_41,Phit_41,~,~,~,~] = FAIR(2,ht_41,f_41,'BE');
Tt_45 = Tt_45R;
k = 0;
k_max = 100;
diff = 0.5;
while diff > 0.01 && k < k_max
[TtT45,PtP45,MFP_45] = MASSFP_NR(Tt_45,f_45,M_45);
pi_tH = sqrt(Tt_45/Tt_4)*A4_A45*m_45*MFP_4/(m_4*MFP_45);
%[~,h,Pr,Phi,~,~,~,~] = FAIR(1,Tt_45R,f_45,'BE'); %??????
Prt_44i = pi_tH*Prt_41;
[~,ht_44i,Pr_44i,Phi_44i,~,~,~,~] = FAIR(3,Prt_44i,f_41,'BE');
ht_44 = ht_41-eta_tH*(ht_41-ht_44i);
tau_tH = ht_44/ht_41;
ht_45 = (m_41*ht_44+e2*ht_3)/(m_41+e2);
[Tt_45n,h,Pr,Phi,~,~,~,~] = FAIR(2,ht_45,f_45,'BE');
diff = abs(Tt_45n-Tt_45);
k = k+1;
Tt_45 = Tt_45n;
end
end
%% MIXER: VSH Mixer Solver
function [M_6A,MFP_6A,Tt_6A,gamma_6A,Cp_6A] = MIXER_VSH(M_6,M_16,gamma_16,gamma_6,R_6,alpha_prime,A16_A6,T_6,f_6A,ht_6A)
C = sqrt(R_6*T_6/gamma_6)*((1+gamma_6*M_6^2)+A16_A6*(1+gamma_16*M_16^2))/(M_6*(1+alpha_prime));
M_p = M_6; % Initially guess that M_6A = M_6.
error = 1;
k = 0;
[Tt_6A,~,Prt_6A,~,~,~,~,~] = FAIR(2,ht_6A,f_6A,'BE');
[TtT6A,~,MFP_6A] = MASSFP_NR(Tt_6A,f_6A,M_p);
while abs(error) > 0.00000001 && k <= 100
    T_6A = Tt_6A/TtT6A;
    [~,h_6A,Pr_6A,Phi_6A,Cp_6A,R_6A,gamma_6A,a_6A] = FAIR(1,T_6A,f_6A,'BE');
    f = sqrt(R_6A*T_6A/gamma_6A)*(1+gamma_6A*M_p^2)/M_p - C; % Non-Linear function to be iterated.
    d = (-1/(M_p^2) + gamma_6A*M_p)*sqrt(R_6A*T_6A/gamma_6A); % Derivative of the Non-Linear function.
    M_6A = M_p - f/d; % New approximation resulting from the iteration.
    M_p = M_6A;
   [TtT6A,~,MFP_6A] = MASSFP_NR(Tt_6A,f_6A,M_p);
    error = abs(M_6A-M_p); % Update the error measured.
    k = k+1;
end
end
%% RGCOMPR:
function [M,TtT,PtP,MFP] = RGCOMPR_NR(mode,Tt,f,IN)
switch mode
    case 1 %Mach Number known
        M = IN;
        [TtT,PtP,MFP] = MASSFP_NR(Tt,f,M);
    case {2,3}
        [~,ht,Prt,Phit,~,~,~,~] = FAIR(1,Tt,f,'BE');
        if mode == 2 %TtT known;
            TtT = IN;
            T = Tt/TtT;
            [~,h,Pr,Phi,Cp,R,gamma,a] = FAIR(1,T,f,'BE');
            PtP = Prt/Pr;
        else %PtP known;
            PtP = IN;
            Pr = Prt/PtP;
            [T,h,~,Phi,Cp,R,gamma,a] = FAIR(3,Pr,f,'BE');
            TtT = Tt/T;
        end
        V2 = 2*31.174*778.16*(ht-h);
        if V2 < 0
            M = 0;
            TtT = 1;
        else
            M = sqrt(V2)/a;
        end
        MFP = (M/PtP)*sqrt(gamma*TtT*32.174/R);
    case {4,5} %MFP Known
        MFP = IN;
        [TtT,PtP,M] = INV_MASSP_NR(mode,Tt,f,MFP);
end
end
%% MIL-5008B: Ram Recovery Model
function eta_R_spec = MIL_E_5008B(M_0)
if M_0 <= 1
    eta_R_spec = 1;
elseif 1 < M_0 && M_0 <5
    eta_R_spec = 1 - 0.075*power(M_0 - 1,1.35);
else
    eta_R_spec = 800/(M_0^4 +935);
end
end