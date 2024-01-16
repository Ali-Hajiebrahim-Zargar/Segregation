clear all
clc;
% pattern:
% substrate-->AlAs
% 6 MLs AlAs 6 MLs AlAs
% growth rate= 0.5 ML/s
% E1=1.7 ev
% E2=1.8 ev
%------  different temperature Ga growth rate=1   ----------
figure
growth_rate=0.5;
ML=40;
ML_counter=1;
discrete=0;
time=-1:0.001:30;
kB=1.38064852*10^-23;
E1=1.7*1.602*10^-19;
E2=1.8*1.602*10^-19;
X_b_Ga_0=0;
X_b_Al_0=1;
X_s_Ga_0=0;
X_s_Al_0=0;
temp=743;
phi_Al=0*growth_rate;
phi_Ga=1*growth_rate;
v=10^13;
p1=v*exp(-E1/(kB*temp));
p2=v*exp(-E2/(kB*temp));
tspan = [0:0.01:1/growth_rate];
y0 = 0;
for c=0:ML+1
ML_counter=ML_counter+1;
Z=p1-p2;
L=-p1*X_s_Al_0-2*p1*X_s_Ga_0-p1*X_b_Ga_0+p2*X_s_Ga_0-p2*X_b_Al_0;
W=-p1*phi_Al-2*p1*phi_Ga+p2*phi_Ga;
G=p1*X_s_Ga_0*phi_Al+2*p1*X_s_Ga_0*phi_Ga+p1*X_s_Al_0*phi_Ga+p1*X_b_Ga_0*phi_Al+p1*X_b_Ga_0*phi_Ga;
T=p1*phi_Ga*phi_Ga+p1*phi_Ga*phi_Al;
M=phi_Ga+p1*X_s_Al_0*X_s_Ga_0+p1*X_s_Ga_0*X_s_Ga_0+p1*X_b_Ga_0*X_s_Al_0+p1*X_b_Ga_0*X_s_Ga_0;
[t,y] = ode45(@(t,y) Z*y*y+L*y+W*y*t+G*t+T*t*t+M, tspan, y0);
XsGa=y(100/growth_rate+1);
XsAl=(phi_Al+phi_Ga)*t(100/growth_rate+1)-XsGa;
XbAl=-XsAl+X_b_Al_0+phi_Al*t(100/growth_rate+1);
XbGa=-XsGa+X_b_Ga_0+phi_Ga*t(100/growth_rate+1);
X_b_Ga_0=XsGa;
X_b_Al_0=XsAl;
X_b_Ga(1,ML_counter)=XbGa;
if (ML_counter>6)&&(ML_counter<13)
phi_Al=1*growth_rate;
phi_Ga=0*growth_rate;  
end
if (ML_counter>12)&&(ML_counter<18)
phi_Al=0*growth_rate;
phi_Ga=1*growth_rate;    
end
discrete=discrete+(X_b_Ga(1,ML_counter)-X_b_Ga(1,ML_counter-1))*heaviside(time-ML_counter+3);
end
ml=[-1:1:ML+1];
[ml;X_b_Ga]
plot(time,discrete,'m-','linewidth',1)
title('GaAs/AlAs growth rate=1 : temp=500 & temp=400')
xlabel('Thickness(ML)')
ylabel('Ga concentration')
axis([-1 30 0 1])
hold on
%-------------------------------------------------------------------------------------------
