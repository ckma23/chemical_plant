close all
clear all
clc

%----------------------------------------------
%This Code is for  for the Design Equations of the PFR Reactor to
%obtain Reactor Concentrations. Assume Isothermal and Isobaric Reactor

global Temp
global Cin_PE
global Cin_Cl2
global MR
%% Notation used for this code
%CH2=CHCH2Cl  : [AC],     "Allyl Chloride"
%CH2ClCHClCH3 : [12DCP],  "1,2 DiChloropropylene"
%CH2=CHCH2Cl  : [13DCP],  "1,3 DiChloropropylene"
%CH2=CHCH3    : [PE],     "Propylene"
%Cl2          : [Cl2],    "Chloride"
%HCl          : [HCl],    "Hydrochloric Acid"

%% Molecular Weight (g/mole)
MW_AC=76.52;
MW_12DCP=112.99;
MW_13DCP=110.97;
MW_PE=42.08;
MW_Cl2=70.91;
MW_HCl=36.4;
% %% Heat of Vaporization
% deltaH13DCP=113; % Btu/lb
% deltaH12DCP=122; % Btu/lb
% deltaHPE=188.28; % Btu/lb
% deltaHAC=159.07; % Btu/lb
% deltaHHCl=178; % Btu/lb
% deltaHCl2=123.73; % Btu/lb
 %% Heat of Vaporization in J/g
 deltaH13DCP=113*2.31; % J/g
 deltaH12DCP=122*2.31; % J/g
 deltaHPE=188*2.31; % J/g
 deltaHAC=159*2.31; % J/g
% deltaHHCl=178; % J/g
% deltaHCl2=123.73; % J/g

 %% Heat of Vaporization in J/mole
 deltaH13DCP=263*MW_13DCP; % J/mole
 deltaH12DCP=283*MW_12DCP; % J/mole
 deltaHPE=437.28*MW_PE; % J/mole
 deltaHAC=369.07*MW_AC; % J/mole
% deltaHHCl=178; % J/mole
% deltaHCl2=123.73; % J/mole


Ptotal=1500*10^3; % Operating Pressure in Pascals

uppertau=.0001; % set upper limit on tau ***optimum***
lowertau=0;  % set lower limit on tau
taus=[lowertau uppertau];

Tems=650; % optimum operating temperature


    % solve the PFR FUN for the PFR Design equations where CTi changes 
for k=1:length(Tems)
  
         Temp=Tems(k);
    
    Ptotal=1500*10^3;                                % operating pressure units of Pascal
    Ptotalatm=Ptotal/(101.3*10^3);                   % operating pressure units of atm
    MR=1;                                            % MR=Cl2/PE
    R=.082;                                          % gas constant in atm*L/mole*K
    Cin_PE=MR/(1+MR);             % inlet concentration of PE in mole/m^3
    Cin_Cl2=1/(1+MR);           % inlet concentration of Cl2 in mole/m^3
    Cin_HCl=0;                   % inlet concentration of HCl
    Cin_12DCP=0;                 % inlet concentration of 12DCP
    Cin_13DCP=0;                 % inlet concentration of 13DCP
    Cin_AC=0;                    % inlet concentration of AC
    
        xin=[Cin_AC;Cin_12DCP;Cin_13DCP;Cin_PE;Cin_Cl2;Cin_HCl];       % inlet concentration vector
         
         [tau,x] = ode45(@PFRFUNCTION,taus,xin,1e-10);
         clc;
    % concentrations are in mole/m^3     
    tempsol_AC(1,:)=x(:,1);                  % store the concentration solutions for ever temperature
    tempsol_12DCP(1,:)=x(:,2);               % store the concentration solutions for ever temperature
    tempsol_13DCP(1,:)=x(:,3);               % store the concentration solutions for ever temperature
    tempsol_PE(1,:)=x(:,4);                  % store the concentration solutions for ever temperature            
    tempsol_Cl2(1,:)=x(:,5);                 % store the concentration solutions for ever temperature 
    tempsol_HCl(1,:)=x(:,6);                 % store the concentration solutions for ever temperature 
    tempsol_tau=tau;

end

pick=1;
% selectivityplot=tempsol_AC(pick,:)./(tempsol_PE(pick,:));
conversionplot_PE=(Cin_PE-tempsol_PE(pick,:))./Cin_PE;
conversionplot_Cl2=(Cin_Cl2-tempsol_Cl2(pick,:))./Cin_Cl2;
%% From level 2 Balance

P_AC=77.79*1000;         % moles/hr
PR=tempsol_13DCP(pick,:)./tempsol_AC(pick,:);
W=tempsol_12DCP(pick,:)./tempsol_AC(pick,:);
P_13DCP=PR.*P_AC;
F_PE=(1+PR+W).*P_AC;
F_Cl2=F_PE+PR.*P_13DCP;
P_HCl=P_AC+2.*P_13DCP;
P_12DCP=W*P_AC;
F_Cl2=F_PE/MR;

% P_13DCP=PR*P_AC;
% F_PE=P_AC./selectivityplot;
% P_HCl=P_AC+2.*P_13DCP;
% P_12DCP=F_PE-P_AC-P_13DCP;
% F_Cl2=P_AC+2.*P_13DCP+P_12DCP;
R=.082;
%% find Rho
rho_PE=Ptotalatm./(R*Temp(pick)).*tempsol_PE(pick,:);
rho_13DCP=Ptotalatm./(R*Temp(pick)).*tempsol_PE(pick,:);
rho_AC=Ptotalatm./(R*Temp(pick)).*tempsol_PE(pick,:);
rho_HCl=Ptotalatm./(R*Temp(pick)).*tempsol_PE(pick,:);
rho_Cl2=Ptotalatm./(R*Temp(pick)).*tempsol_PE(pick,:);
rho_12DCP=Ptotalatm./(R*Temp(pick)).*tempsol_PE(pick,:);

v0=(F_Cl2./rho_Cl2)+(F_PE./rho_PE); %moles/s divided by mole/L
V=tau.*v0'; %in liters



figure(1)
hold on
plot(tau,tempsol_AC(1,:),'red')
hold on
plot(tau,tempsol_12DCP(1,:),'blue')
hold on
plot(tau,tempsol_13DCP(1,:),'green')
hold on
plot(tau,tempsol_PE(1,:),'black')
hold on
plot(tau,tempsol_Cl2(1,:),'cyan')
hold on
plot(tau,tempsol_HCl(1,:),'magenta')
xlabel('Volume'),ylabel('Dimensionless Concentration'),legend('AC','12DCP','13DCP','PE','Cl2','HCl')
 
% figure(2)
%  hold on
%  plot(conversionplot,V,'red')
%  
 %% Fresh Feed Streams
% Flow of PE into the reactor
%F_PE(end)
% Flow of PE out the reactor is
 Fout_PE=((1-conversionplot_PE(end))*F_PE(end));
% Flow of Cl2 into the reactor
%F_Cl2(end)
% Flow of Cl2 out the reactor is
 Fout_Cl2=(1-conversionplot_Cl2(end))*F_Cl2(end);
 FFin_PE=F_PE(end)-Fout_PE;
 FFin_Cl2=F_Cl2(end)-Fout_Cl2;

  fprintf('The Temperature of the reactor is %g K.\n',Temp(pick))
  fprintf('The Pressure of the system is %g Pascals.\n',Ptotal)
  fprintf('The Molar Ratio of the system is %g .\n',MR)
  fprintf('The conversion at the end of the reactor for Cl2 is %g .\n', conversionplot_Cl2(end))
  fprintf('The conversion at the end of the reactor for PE is %g .\n', conversionplot_PE(end))
%   fprintf('The selectivity at the end of the reactor is %g .\n', selectivityplot(end))
  fprintf('The flow of PE into the reactor is %g mole/hr or %g g/hr.\n',F_PE(end),F_PE(end)*MW_PE)
  fprintf('The flow of Cl2 into the reactor is % g mole/hr or %g g/hr.\n',F_Cl2(end),F_Cl2(end)*MW_Cl2)
   fprintf('The flow of PE out the reactor is %g mole/hr or %g g/hr to be recycled .\n',Fout_PE,Fout_PE*MW_PE)
   fprintf('The flow of Cl2 out the reactor is %g mole/hr or %g g/hr to be recycled .\n',Fout_Cl2,Fout_Cl2*MW_Cl2)
  fprintf('The flow of AC out of the reactor is %g mole/hr or %g g/hr  .\n',P_AC(end),P_AC(end)*MW_AC)
  fprintf('The flow of 12DCP out of the reactor is %g mole/hr or %g g/hr .\n',P_12DCP(end),P_12DCP(end)*MW_12DCP)
  fprintf('The flow of 13DCP out of the reactor is %g mole/hr or %g g/hr .\n',P_13DCP(end),P_13DCP(end)*MW_13DCP)
  fprintf('The flow of HCl out of the reactor is %g mole/hr or %g g/hr .\n',P_HCl(end),P_HCl(end)*MW_HCl)
   fprintf('The flow of fresh PE is %g mole/hr  or %g g/hr.\n',FFin_PE,FFin_PE*MW_PE)
   fprintf('The flow of fresh Cl2 is %g mole/hr  or %g g/hr.\n',FFin_Cl2,FFin_Cl2*MW_Cl2)
  fprintf('The flow "vo" required is %g L/hr.\n',v0(end)*3600)
  fprintf('The size of the PFR reactor is %g L or %g m^3.\n',V(end),V(end)/1000)
  fprintf('The product ratio of 13DCP to AC is %g.\n',PR(end)) 
  fprintf('The waste ratio of 12DCP to AC is %g .\n',W(end))
  %% Assume HCl and Cl2 is removed before Botillation
  fprintf('The flow rate "Feed" into the separation system is %g mole/hr .\n',P_12DCP(end)+P_13DCP(end)+P_AC(end)+Fout_PE(end))
  fprintf('The mole fraction of PE in the Feed into the separation system is %g .\n',Fout_PE(end)/(P_12DCP(end)+P_13DCP(end)+P_AC(end)+Fout_PE(end)))
  fprintf('The mole fraction of AC in the Feed into the separation system is %g .\n',P_AC(end)/(P_12DCP(end)+P_13DCP(end)+P_AC(end)+Fout_PE(end)))
 fprintf('The mole fraction of 12DCP in the Feed into the separation system is %g .\n',P_12DCP(end)/(P_12DCP(end)+P_13DCP(end)+P_AC(end)+Fout_PE(end)))
  fprintf('The mole fraction of 13DCP in the Feed into the separation system is %g .\n',P_13DCP(end)/(P_12DCP(end)+P_13DCP(end)+P_AC(end)+Fout_PE(end)))
%  fprintf('The mole fraction of HCl in the Feed into the separation system is %g .\n',(tempsol_HCl(pick,end)*v0)/((tempsol_PE(pick,end)+tempsol_AC(pick,end)+tempsol_12DCP(pick,end)+tempsol_13DCP(pick,end)+tempsol_HCl(pick,end)+tempsol_Cl2(pick,end))*v0))
%  fprintf('The mole fraction of Cl2 in the Feed into the separation system is %g .\n',(tempsol_Cl2(pick,end)*v0)/((tempsol_PE(pick,end)+tempsol_AC(pick,end)+tempsol_12DCP(pick,end)+tempsol_13DCP(pick,end)+tempsol_HCl(pick,end)+tempsol_Cl2(pick,end))*v0))





%% Costing Section

% Raw material costs
CPE=1.87; % $/kg -- cost of propylene
CCl2=0.36; % $/kg -- cost of chlorine
CHCl=0.085; % $/kg -- cost of HCl
CAC=2.53; % $/kg -- cost of allyl chloride
CFG=10; % $/gal -- cost of Floragrow
C12DCP=0.2; % $/kg, must incinerate -- 1,2 Dichloropropane
% incineration cost is between 20-40 cents: chosen to be 20 cents

%% Finance Parameters
Ncon=2; % Nconstruction
Nop=10; % Noperations
ratetax=0.48; % Tax rate percent
ratecon=0.06; % Construction rate percent
ratefin=0.04; % Finance rate percent
rateent=0.08; % Enterprise rate percent
a=[0.5 0.5 0 0];
b=[0.8 0.9 0.95 1 1 1 1 1 1 1];
alphaWC=0.2; % alpha value for working capital
alphaSU=0.1; % alpha value for start up capital
alphaSV=0.03; % alpha value for salvage value

%% Costing Parameters
MS=1600; % assumed M&S index 2013
laborcost=1800000; % $/yr -- labor cost

%% Conditions for the reactor
Qh=(1.327*10^7)*(1000/1055); 
QR=-35190; % kJ/hr -- heat flow into the reactor
Vfeet=V(end)*(1/28.3168466); % volume of reactor in ft^3
Dfeet=((4*Vfeet)/(3.14159*6))^(1/3); % diameter of reactor in ft^3
Lfeet=(4*Vfeet)/(3.14159*(Dfeet^2)); % length of reactor in ft^3
number_of_reactors=1;

% Stripping column
%
%% OBTAIN VALUES FROM ASPEN PLUS
Feedtosep=P_12DCP(end)+P_13DCP(end)+P_AC(end)+Fout_PE(end); % g mole/hr

%% FIRST Separation System PULLING OF PE
% fprintf('The Distillate Flowrate out of Sep 1 of Benzene is %g mole/hr .\n',temp_solB(pick,end)*v0)
% fprintf('The Bottom Flowrate out of Sep 1 of Benzene is %g mole/hr .\n',temp_solT(pick,end)*v0+temp_solP(pick,end)*v0+temp_solTMB(pick,end)*v0)
% fprintf('The xA out of Sep 1 is %g mole/hr .\n',(temp_solB(pick,end)*v0)/(temp_solB(pick,end)*v0))
% fprintf('The xD out of Sep 1 is %g mole/hr .\n',(temp_solTMB(pick,end)*v0)/(temp_solT(pick,end)*v0+temp_solP(pick,end)*v0+temp_solTMB(pick,end)*v0))
q=1;
% T=; operating temperature of the column
% P= 15; bar,operating pressure of the column 

Fsep1=Feedtosep; %Feed into column 1 %g mole/hr
Dsep1=0.02948545*1000*3600; %Distillate from column 1 % gmole/hr
Bsep1=0.03293316*1000*3600; %Bottoms from column 1 % gmole/hr
Rsep1=.6; % reflux ratio
Ssep1=2.37; % reboil ratio

Nrealsep1=14.4; % number of columns 
Dist1_PEflow=0.02945596*1000*3600;% gmole/hr
Dist1_ACflow=2.1279e-05*1000*3600;% gmole/hr
Dist1_12DCPflow=8.2036e-06*1000*3600;% gmole/hr
Dist1_13DCPflow=2.9485e-09*1000*3600;% gmole/hr
XDist1_PE=Dist1_PEflow/(Dist1_PEflow+Dist1_12DCPflow+Dist1_13DCPflow);
XDist1_AC=Dist1_ACflow/(Dist1_PEflow+Dist1_12DCPflow+Dist1_13DCPflow);
XDist1_12DCP=Dist1_12DCPflow/(Dist1_PEflow+Dist1_12DCPflow+Dist1_13DCPflow);
XDist1_13DCP=Dist1_13DCPflow/(Dist1_PEflow+Dist1_12DCPflow+Dist1_13DCPflow);
Bot1_PEflow=3.2933e-10*1000*3600;% gmole/hr
Bot1_ACflow=0.02158710*1000*3600;% gmole/hr
Bot1_12DCPflow=0.00832243*1000*3600;% gmole/hr
Bot1_13DCPflow=0.00302361*1000*3600;% gmole/hr
XBot1_PE=Bot1_PEflow/(Bot1_PEflow+Bot1_12DCPflow+Bot1_13DCPflow);
XBot1_AC=Bot1_ACflow/(Bot1_PEflow+Bot1_12DCPflow+Bot1_13DCPflow);
XBot1_12DCP=Bot1_12DCPflow/(Bot1_PEflow+Bot1_12DCPflow+Bot1_13DCPflow);
XBot1_13DCP=Bot1_13DCPflow/(Bot1_PEflow+Bot1_12DCPflow+Bot1_13DCPflow);

VBsep1=Ssep1*Bsep1; % mol/hr
VTsep1=(Rsep1+1)*Dsep1; % mol/hr
lambda1D=XDist1_PE*deltaHPE+XDist1_13DCP*deltaH13DCP+XDist1_12DCP*deltaH12DCP+XDist1_AC*deltaHAC; %J/mole
lambda1B=XBot1_PE*deltaHPE+XBot1_13DCP*deltaH13DCP+XBot1_12DCP*deltaH12DCP+XBot1_AC*deltaHAC; %J/mole
enthalpyofdist1=7558463.72;%J/kmole
enthalpyofbot1=-29586656;%J/kmole
QCsep1=lambda1D*VTsep1; %  Joules/hr
QRsep1=lambda1B*VBsep1; % Joules/hr
Fcoolsep1=QCsep1/(lambda1D*15); % mol/hr
  

%% Second SEPARATION SYSTEM PULLING OFF 13DCP
Fsep2=Bsep1; %Feed into column 2
Dsep2=0.03015554*1000*3600; %Distillate out of column 2
Bsep2=0.00277761*1000*3600; %Bottoms out of column 2
Rsep2=8.5; % reflux ratio
Ssep2=103.01; % reboil ratio;
Nrealsep2=50.9; % number of stages

Dist2_PEflow=0*1000*3600;% gmole/hr
Dist2_ACflow=0.02158710*1000*3600;% gmole/hr
Dist2_12DCPflow=0.00826688*1000*3600;% gmole/hr
Dist2_13DCPflow=0.00030155*1000*3600;% gmole/hr
XDist2_PE=Dist2_PEflow/(Dist2_PEflow+Dist2_12DCPflow+Dist2_13DCPflow);
XDist2_AC=Dist2_ACflow/(Dist2_PEflow+Dist2_12DCPflow+Dist2_13DCPflow);
XDist2_12DCP=Dist2_12DCPflow/(Dist2_PEflow+Dist2_12DCPflow+Dist2_13DCPflow);
XDist2_13DCP=Dist2_13DCPflow/(Dist2_PEflow+Dist2_12DCPflow+Dist2_13DCPflow);
Bot2_PEflow=0*1000*3600;% gmole/hr
Bot2_ACflow=2.7776e-14*1000*3600;% gmole/hr
Bot2_12DCPflow=5.5552e-05*1000*3600;% gmole/hr
Bot2_13DCPflow=0.00272206*1000*3600;% gmole/hr
XBot2_PE=Bot2_PEflow/(Bot2_PEflow+Bot2_12DCPflow+Bot2_13DCPflow);
XBot2_AC=Bot2_ACflow/(Bot2_PEflow+Bot2_12DCPflow+Bot2_13DCPflow);
XBot2_12DCP=Bot2_12DCPflow/(Bot2_PEflow+Bot2_12DCPflow+Bot2_13DCPflow);
XBot2_13DCP=Bot2_13DCPflow/(Bot2_PEflow+Bot2_12DCPflow+Bot2_13DCPflow);

 
 VBsep2=Ssep2*Bsep2; %mol/hr
 VTsep2=(Rsep2+1)*Dsep2; %mol/hr

 lambda2D=XDist2_PE*deltaHPE+XDist2_13DCP*deltaH13DCP+XDist2_12DCP*deltaH12DCP+XDist2_AC*deltaHAC; %J/mole
 lambda2B=XBot2_PE*deltaHPE+XBot2_13DCP*deltaH13DCP+XBot2_12DCP*deltaH12DCP+XBot2_AC*deltaHAC; %J/mole
 enthalpyofdist2=-9497823.9;%J/kmole
 enthalpyofbot2=-64819174;%J/kmole
 QCsep2=lambda2D*VTsep2; %Joule/hr
 QRsep2=lambda2B*VBsep2; %Joule/hr
 Fcoolsep2=QCsep2/(lambda2D*15); % mol/hr

%% Third SEPARATION SYSTEM AC to 12DCP
Fsep3=(Dist2_ACflow+Bot2_ACflow+Dist2_12DCPflow+Bot2_12DCPflow)*1000*3600;%Feed into column 3
Dsep3=Dist2_ACflow+Bot2_ACflow*1000*3600;%Distillate column 3
Bsep3=(Dist2_12DCPflow+Bot2_12DCPflow)*1000*3600 ;%Bottoms column 3
alphaACto12DCP=2.6;

xsep3=((Dist2_ACflow+Bot2_ACflow)*1000*3600)/Fsep3; % mole fraction of AC to 12DCP
rminsep3=1/((alphaACto12DCP-1)*xsep3); %minimum reflux
Rsep3=1.5*rminsep3; % reflux ratio
Ssep3=(Dsep3/Bsep3)*(Rsep3+q)-(1-q);
xbsep3=.003;
xdsep3=.997;
Nminsep3=log((xdsep3/(1-xdsep3)*((1-xbsep3)/xbsep3))/log(alphaACto12DCP));
Nsep3=(.75*(1-((Rsep3-rminsep3)/(Rsep3+1))^.5688)+Nminsep3)/(1-.75*(1-((Rsep3-rminsep3)/(Rsep3+1))^.5688));
a=.24;
mu0=1*10^-3;
musep3=1*10^-4;
E0sep3=exp(-sqrt(alphaACto12DCP*(musep3/mu0)))*(1-a)+a;
Nrealsep3=Nsep3/E0sep3;
VBsep3=Ssep3*Bsep3; %mol/hr
VTsep3=(Rsep3+1)*Dsep3; %mol/hr
lambda3D=deltaHAC;
lambda3B=deltaH12DCP; %joule/mol
QCsep3=lambda3D*VTsep3; % in joule/hr
QRsep3=lambda3B*VBsep3; % in joule/hr
Fcoolsep3=QCsep3/(lambda3D*15); % mol/hr
Fcool=(Fcoolsep1+Fcoolsep2+Fcoolsep3)*.018; %kg/hr

% Separation terms
Ht=0.6096; % Spacing between trays, m
heiA=(3*Ht)+(Ht*Nrealsep1);
heiB=(3*Ht)+(Ht*Nrealsep2);
heiC=(3*Ht)+(Ht*Nrealsep3);
Hca=heiA*3.28; % Convert column A height to feet
Hcb=heiB*3.28; % Convert column B height to feet
Hcc=heiC*3.28; % Convert column C height to feet

MaveABCDEF=(MW_AC+MW_12DCP+MW_13DCP+MW_PE+MW_Cl2+MW_HCl)/6; % average molecular weight of the 6 compounds
rholPE=613.9; % Liquid density of propylene, kg/m^3
rholCl2=1562.5; % Liquid density of chlorine, kg/m^3
rholHCl=1191; % Liquid density of hydrogen chloride, kg/m^3
rholAC=940; % Liquid density of allyl chloride, kg/m^3
rhol13DCP=1220; % Liquid density (avg) of trans & cis-1,3-dichloropropene, kg/m^3
rhol12DCP=1156; % Liquid density of 1,2-dichloropropane, kg/m^3
% MaveABCD=(BMM+PMM+TMM+TMBMM)/4 ; % average molecular weight of benzene, toluene, paraxylene, and TMB
% MaveBCD=(PMM+TMM+TMBMM)/3; % average molecular weight of toluene, paraxylene, and TMB
% MaveCD=(PMM+TMBMM)/2; % average molecular weight of paraxylene and TMB
% rholB=873.8; % liquid density kg/m^3
% rholP=861; % liquid density kg/m^3
% rholT=867; % liquid density kg/m^3
% rholTMB=876.1; % liquid density kg/m^3
Tboil13DCP=103.8+273; % boiling point 13 DCP, K
Tboil12DCP=96.8+273; % boiling point 12DCP, K
TboilPE=-47+273; % boiling point PE, K
TboilCl2=-33.7+273; % boiling point Cl2, K
TboilHCl=-83.7+273; % boiling point HCl, K
TboilAC=44.9+273; % boiling point AC, K

% rhovA_top=(10*BMM)/(0.0821*TboilB); % liquid density, column A top, kg/m^3
% rhovA_bottom=(10*((PMM+TMM+TMBMM)/3))/(0.0821*TboilT); % liquid density, column A bottom, kg/m^3
% rhovB_top=(10*TMM)/(0.0821*TboilT); % liquid density, column B top, kg/m^3
% rhovB_bottom=(10*((PMM+TMBMM)/2))/(0.0821*TboilP); % liquid density, column B bottom, kg/m^3
% rhovC_top=(10*PMM)/(0.0821*TboilP); % liquid density, column C top, kg/m^3
% rhovC_bottom=(10*TMBMM)/(0.0821*TboilTMB); % liquid density, column C bottom, kg/m^3
% VA=(VBsep1+VTsep1)/2;
% VB=(VBsep2+VTsep2)/2;
% VC=(VBsep3+VTsep3)/2;
% diaA=2*sqrt((MaveABCD*1.25*VA)/(3.14159*sqrt(rhovA_top*rhovA_bottom)*0.6*439));
% diaB=2*sqrt((MaveBCD*1.25*VB)/(3.14159*sqrt(rhovB_top*rhovB_bottom)*0.6*439));
% diaC=2*sqrt((MaveCD*1.25*VC)/(3.14159*sqrt(rhovC_top*rhovC_bottom)*0.6*439));
% Dca=diaA*3.28;
% Dcb=diaB*3.28;
% Dcc=diaC*3.28;
%     ARBC=237;
%     ABC=571;
%     xBsep1T=((temp_solT(pick,end)*v0)/(temp_solT(pick,end)*v0+temp_solP(pick,end)*v0+temp_solTMB(pick,end)*v0));
%     xBsep1P=((temp_solP(pick,end)*v0)/(temp_solT(pick,end)*v0+temp_solP(pick,end)*v0+temp_solTMB(pick,end)*v0));
%     xBsep1TMB=((temp_solTMB(pick,end)*v0)/(temp_solT(pick,end)*v0+temp_solP(pick,end)*v0+temp_solTMB(pick,end)*v0));
%     xBsep2P=v0*temp_solP(pick,end)/(temp_solP(pick,end)*v0+temp_solTMB(pick,end)*v0);
%     xBsep2TMB=v0*temp_solTMB(pick,end)/(temp_solP(pick,end)*v0+temp_solTMB(pick,end)*v0);
%     xBsep3TMB=(v0*temp_solTMB(pick,end))/(temp_solTMB(pick,end)*v0);
%     xCsep1B=(v0*temp_solB(pick,end))/(v0*temp_solB(pick,end));
%     xCsep2T=(v0*temp_solT(pick,end))/(v0*temp_solT(pick,end));
%     xCsep3P=(v0*temp_solP(pick,end))/(v0*temp_solP(pick,end));
%     AR1=(1/11250)*.000947*(VTsep1*1000*(xBsep1T*deltaHT+xBsep1P*deltaHP+xBsep1TMB*deltaHTMB)); %square ft
%     AR2=(1/11250)*.000947*(VBsep2*1000*(xBsep2P*deltaHP+xBsep2TMB*deltaHTMB));%square ft
%     AR3=(1/11250)*.000947*(VBsep3*1000*(xBsep3TMB*deltaHTMB)); %square ft
%     TboilB=(180.7*(9/5))+32; %in Farenheit
%     TboilT=(217.6*(9/5))+32; %in Farenheit
%     TboilP=(253.1*(9/5))+32; %in Farenheit
%     TboilTMB=(287.7*(9/5))+32; %in Farenheit
%     AC1=VTsep1*(xCsep1B)*.000947*deltaHB*1000*(1/3000)*(log((TboilB-90)/(TboilB-120)));
%     AC2=VBsep2*(xCsep2T)*.000947*deltaHT*1000*(1/3000)*(log((TboilT-90)/(TboilT-120)));
%     AC3=VBsep3*(xCsep3P)*.000947*deltaHP*1000*(1/3000)*(log((TboilP-90)/(TboilP-120)));
% %     Costreb=(MS/280)*(101.3*10)*(237^.65)*(3.29/3);
% %     Costcond=(MS/280)*(101.3*10)*(571^.65)*(3.29/3);
%     Costreb=(MS/280)*(101.3)*(237^.65)*(3.29/3);
%     Costcond=(MS/280)*(101.3)*(571^.65)*(3.29/3);

%% Separation System Constants from Douglas 
% Separation constants
Fmc=1;
Fpc=1.2;
Fst=1;
Ftt=0;
Fmt=0;
Fct=Fst+Ftt+Fmt;
Fmh=1;
Fdh=0.8;
Fph=0.1;
Fch=(Fdh+Fph)*Fmh;
CWcost=0.00008; % cost of cooling water, $/kg
a=1:5
 for i=1:length(Temp)
    

    % Heater cost
    Value(1,i)=(MS/280)*5.07*10^3*((Qh/(10^6))^0.85)*1; % Purchase cost, $
    Value(2,i)=(MS/280)*5.07*10^3*((Qh/(10^6))^0.85)*2.23; % Install cost, $
    
    % Reactor cost
    Value(3,i)=(MS/280)*101.9*(Dfeet^1.066)*(Lfeet^0.82)*((2.18+(1*1.2))/3)*number_of_reactors; % Annual cost, $/yr
    Value(4,i)=(MS/280)*101.9*(Dfeet^1.066)*(Lfeet^0.82)*(1*1.2)*number_of_reactors; % Purchase cost, $
    Value(5,i)=(MS/280)*101.9*(Dfeet^1.066)*(Lfeet^0.82)*(2.18+(1*1.2))*number_of_reactors; % Install cost, $
    
    % Separator costs (Trays and Shells)
%     Value(6,i,ii)=(MS/280)*(101.9)*(Dca^1.066)*(Hca^0.82)*(Fmc*Fpc); % Purchase cost of shell (Column A)
%     Value(7,i,ii)=(MS/280)*(101.9)*(Dcb^1.066)*(Hcb^0.82)*(Fmc*Fpc); % Purchase cost of shell (Column B)
%     Value(8,i,ii)=(MS/280)*(101.9)*(Dcc^1.066)*(Hcc^0.82)*(Fmc*Fpc); % Purchase cost of shell (Column C)
%     Value(9,i,ii)=(MS/280)*(101.9)*(Dca^1.066)*(Hca^0.82)*(2.18+(Fmc*Fpc)); % Install cost of shell (Column A)
%     Value(10,i,ii)=(MS/280)*(101.9)*(Dcb^1.066)*(Hcb^0.82)*(2.18+(Fmc*Fpc)); % Install cost of shell (Column B)
%     Value(11,i,ii)=(MS/280)*(101.9)*(Dcc^1.066)*(Hcc^0.82)*(2.18+(Fmc*Fpc)); % Install cost of shell (Column C)
%     Value(12,i,ii)=(MS/280)*(4.7)*(Dca^1.55)*(Hca)*(Fct); % Install tray Separator (Column A), $
%     Value(13,i,ii)=(MS/280)*(4.7)*(Dcb^1.55)*(Hcb)*(Fct); % Install tray Separator (Column B), $
%     Value(14,i,ii)=(MS/280)*(4.7)*(Dcc^1.55)*(Hcc)*(Fct); % Install tray Separator (Column C), $
%
% Must find the following variables before tray and shell costs can be
% done: Dca, Dcb, Dcc, Hca, Hcb, Hcc

    % Separator costs (Reboilers and Condensors)
%     Value(15,i,ii)=((MS/280)*(101.9)*(Aca^0.65)*(3.29/3))+((MS/280)*(101.9)*(Ara^0.65)*(3.29/3)); % Annual cost of condensor and reboiler (Column A), $/yr
%     Value(16,i,ii)=((MS/280)*(101.9)*(Acb^0.65)*(3.29/3))+((MS/280)*(101.9)*(Arb^0.65)*(3.29/3)); % Annual cost of condensor and reboiler (Column B), $/yr
%     Value(17,i,ii)=((MS/280)*(101.9)*(Acc^0.65)*(3.29/3))+((MS/280)*(101.9)*(Arc^0.65)*(3.29/3)); % Annual cost of condensor and reboiler (Column C), $/yr
%
% Must find the following variables before reboiler and condensor annual
% costs can be done: Aca, Acb, Acc, Ara, Arb, Arc.  See douglas if you need help on finding the
% other variables to solve for these.
    
    % Separator costs (Heat Exchangers)
%     Value(18,i,ii)=((MS/280)*(101.3)*(Ara^0.65)*(2.29+Fch));% install cost of Heat Exchanger on reboiler sep A
%     Value(19,i,ii)=((MS/280)*(101.3)*(Aca^0.65)*(2.29+Fch));% install cost of Heat Exchanger on condensor sep A
%     Value(20,i,ii)=((MS/280)*(101.3)*(Arb^0.65)*(2.29+Fch));% install cost of Heat Exchanger on reboiler sep B
%     Value(21,i,ii)=((MS/280)*(101.3)*(Acb^0.65)*(2.29+Fch));% install cost of Heat Exchanger on condensor sep B
%     Value(22,i,ii)=((MS/280)*(101.3)*(Arc^0.65)*(2.29+Fch));% install cost of Heat Exchanger on reboiler sep C
%     Value(23,i,ii)=((MS/280)*(101.3)*(Acc^0.65)*(2.29+Fch));% install cost of Heat Exchanger on condensor sep C
%     Value(24,i,ii)=((MS/280)*(101.3)*(Ara^0.65)*(Fch));% purchase cost of Heat Exchanger on reboiler sep A
%     Value(25,i,ii)=((MS/280)*(101.3)*(Aca^0.65)*(Fch));% purchase cost of Heat Exchanger on condensor sep A
%     Value(26,i,ii)=((MS/280)*(101.3)*(Arb^0.65)*(Fch));% purchase cost of Heat Exchanger on reboiler sep B
%     Value(27,i,ii)=((MS/280)*(101.3)*(Acb^0.65)*(Fch));% purchase cost of Heat Exchanger on condensor sep B
%     Value(28,i,ii)=((MS/280)*(101.3)*(Arc^0.65)*(Fch));% purchase cost of Heat Exchanger on reboiler sep C
%     Value(29,i,ii)=((MS/280)*(101.3)*(Acc^0.65)*(Fch));% purchase cost of Heat Exchanger on condensor sep C
%
% See annual costs above.

    % Catalyst cost
    Value(30,i)=0; % cost, $
    
    % Steam cost (due to reboilers)
    Value(31,i)=0; % $/yr
    
    % Cooling Water cost (due to condensors)
    Value(32,i)=CWcost; % $/yr
    
    % Heat Exchanger cost
    Value(33,i)=0; % Purchase cost, $
    
    pcbe(i)=Value(1,i)+Value(4,i)+Value(6,i)+Value(7,i)+Value(8,i)+Value(24,i)+Value(25,i)+Value(26,i)+Value(27,i)+Value(28,i)+Value(29,i)+Value(30,i)+Value(33,i); % PCBE cost base equipment per yr
    isbl(i)=Value(2,i)+Value(5,i)+Value(9,i)+Value(10,i)+Value(11,i)+Value(12,i)+Value(13,i)+Value(14,i)+Value(18,i)+Value(19,i)+Value(20,i)+...
        Value(21,i)+Value(22,i)+Value(23,i); % ISBL cost (installed equipment costs)
    R1(i)=(P_AC(end)*MW_AC*(1/1000)*8400*CAC)+(P_13DCP(end)*MW_13DCP*8400*(1/1.217)*(1/3785.4)*CFG)+(P_HCl(end)*MW_HCl*8400*(1/1000)*CHCl); % $/yr revenue
    C1(i)=((FFin_Cl2*MW_Cl2)*(1/1000)*8400*CCl2)+((FFin_PE*MW_PE)*(1/1000)*8400*CPE)+(P_12DCP(end)*MW_12DCP*(1/1000)*8400*C12DCP); % $/yr cost
    PBT(i)=R1(i)-C1(i)-Value(3,i)-Value(15,i)-Value(16,i)-Value(17,i)-Value(31,i)-Value(32,i)-laborcost; % earnings in $/yr
    WCestimate(i)=(((FFin_Cl2*MW_Cl2)*(1/1000)*8400*CCl2)+((FFin_PE*MW_PE)*(1/1000)*8400*CPE))*(2/12); % $/yr -- estimated working capital (2 months of raw material costs)
    FC(i)=(2.28*(isbl(i))); % $ -- fixed capital 
    TI(i)=FC(i)*(1+alphaSU+alphaWC); % $ -- total investment
    
    SU=-FC(i)*alphaSU; % Start-up capital
    WC=-FC(i)*alphaWC; % Working capital
    
    % this loop is designed to calculate the pre-operations numbers in
    % reverse order
    for k=1:4
        PreDF(k)=(1+ratecon)^(k-1); % Discount factors in years 0 to -3
        PreFC(k)=-a(k)*FC(i); % Fixed capital in years 0 to -3
        PreDCF(k)=PreDF(k)*PreFC(k); %Discounted Cash Flows in years 0 to -3
    end
    TCO=PreDCF(1)+PreDCF(2)+PreDCF(3)+PreDCF(4)+SU+WC;
    TCI=-TCO; % Total Capital Investment
    
    for j=1:10 %in years 1 through 10
        BF(j)=-ratefin*TCI; % Bond Financing in year j
        earnings(j)=b(j).*PBT(i); % earnings in year j
        Dep(j)=-0.1.*(1+alphaSU).*FC(i); % depreciation in year j
        Profit(j)=earnings(j)+BF(j);
        PATyearly(j)=(1-ratetax).*(earnings(j)+BF(j)+Dep(j)); % profit after taxes in year j
        CF(j)=PATyearly(j)-Dep(j); % cash flow in year j
        DF(j)=(1+rateent)^(-j); % discount factor in year j
        DCF(j)=CF(j).*DF(j); % discounted cash flow in year j
    end
    
    DCFWCten=DF(10)*FC(i)*alphaWC; % Discounted cash flow of WC year 10
    DCFSVten=DF(10)*FC(i)*alphaSV*(1-ratetax); % Discounted cash flow of SV year 10
    DCFPOTCIten=-TCI*DF(10); % Discounteed cash flow of Pay-off TCI
    NPV0=DCF(1)+DCF(2)+DCF(3)+DCF(4)+DCF(5)+DCF(6)+DCF(7)+DCF(8)+DCF(9)+DCF(10)+...
        DCFWCten+DCFSVten+DCFPOTCIten;
    NPVproj=NPV0/((1+rateent)^Ncon); % NPVproj
    NPV0avg=NPV0/(TCI*Nop);
    NPVprojavg=NPVproj/(TCI*(Nop+Ncon));
    ROI_BT=earnings(10)/TI(i); % ROI_BT percent return on investment
    POT=(1.1.*FC(i))/(((1-ratetax)*sum(PBT(i)))-(ratetax*Dep(j))); % Pay-Out Time (yr)
    
end

PATyearly(end)