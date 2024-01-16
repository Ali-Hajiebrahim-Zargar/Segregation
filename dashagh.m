clear all
clc;
% pattern:
% substrate-->GaAs
% 10 MLs GaAs 2MLs InAs
% growth rate= 1 ML/s
% growth temperature = 400 or 500 c
% E1=1.8 ev
% E2=2 ev
%------  different temperature in growth rate=1   ----------
figure
growth_rate=1;
ML=12;
ML_counter=3;
discrete=0;
time=-1:0.001:13;
kB=1.38064852*10^-23;
E1=1.8*1.602*10^-19;
E2=2*1.602*10^-19;
X_b_In_0=0;
X_b_Ga_0=1;
X_s_In_0=0;
X_s_Ga_0=0;
temp=773.15;
phi_Ga=0*growth_rate;
phi_In=1*growth_rate;
v=10^13;
p1=v*exp(-E1/(kB*temp));
p2=v*exp(-E2/(kB*temp));
tspan = [0:0.01:1/growth_rate];
y0 = 0;
for c=0:ML+1
ML_counter=ML_counter+1;
Z=p1-p2;
L=-p1*X_s_Ga_0-2*p1*X_s_In_0-p1*X_b_In_0+p2*X_s_In_0-p2*X_b_Ga_0;
W=-p1*phi_Ga-2*p1*phi_In+p2*phi_In;
G=p1*X_s_In_0*phi_Ga+2*p1*X_s_In_0*phi_In+p1*X_s_Ga_0*phi_In+p1*X_b_In_0*phi_Ga+p1*X_b_In_0*phi_In;
T=p1*phi_In*phi_In+p1*phi_In*phi_Ga;
M=phi_In+p1*X_s_Ga_0*X_s_In_0+p1*X_s_In_0*X_s_In_0+p1*X_b_In_0*X_s_Ga_0+p1*X_b_In_0*X_s_In_0;
[t,y] = ode45(@(t,y) Z*y*y+L*y+W*y*t+G*t+T*t*t+M, tspan, y0);
XsIn=y(100/growth_rate+1);
XsGa=(phi_Ga+phi_In)*t(100/growth_rate+1)-XsIn;
XbGa=-XsGa+X_b_Ga_0+phi_Ga*t(100/growth_rate+1);
XbIn=-XsIn+X_b_In_0+phi_In*t(100/growth_rate+1);
X_b_In_0=XsIn;
X_b_Ga_0=XsGa;
X_b_In(1,ML_counter)=XbIn;
if ML_counter>4
phi_Ga=1;
phi_In=0;    
end
discrete=discrete+(X_b_In(1,ML_counter)-X_b_In(1,ML_counter-1))*heaviside(time-ML_counter+3);
end
ml=[-1:1:ML+1];
plot(time,discrete,'m-','linewidth',1)
title('InAs/GaAs growth rate=1 : temp=500 & temp=400')
xlabel('Thickness(ML)')
ylabel('In concentration')
axis([-1 13 0 1.5])
hold on
R=0.4;
x0=1;
ML=60;
time=-20:0.001:40;
ML_counter=1;
X(1,1)=0;
for c=0:ML+1
    if ML_counter<3
    X(1,ML_counter+1)=(x0)*(1-R^(ML_counter));  
    end
    if 2<ML_counter
    X(1,ML_counter+1)=(x0)*(1-R^2)*R^(ML_counter-2);  
    end
    ML_counter=ML_counter+1;
end
ml=[-1:1:ML+1];
plot(X,'b-','linewidth',1)
title('2ML InAs/10ML GaAs:gr=1ML/s ,gt=500' )
xlabel('Thickness(ML)')
ylabel('In concentration')
axis([0 15 -0.05 1])
hold on
R=0.5;
x0=1;
ML=60;
time=-20:0.001:40;
ML_counter=1;
X(1,1)=0;
for c=0:ML+1
    if ML_counter<3
    X(1,ML_counter+1)=(x0)*(1-R^(ML_counter));  
    end
    if 2<ML_counter
    X(1,ML_counter+1)=(x0)*(1-R^2)*R^(ML_counter-2);  
    end
    ML_counter=ML_counter+1;
end
ml=[-1:1:ML+1];
plot(X,'k-','linewidth',1)
legend('Kinetic:E1=1.8 , E2=2','Muraki:R=0.4','Muraki:R=0.5')