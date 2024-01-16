clear all
clc;
%------  different temperature in growth rate=1   ----------
growth_rate=1;
ML=20;
ML_counter=3;
discrete=0;
time=-1:0.001:20;
kB=1.38064852*10^-23;
E1=1.6*1.602*10^-19;
E2=1.75*1.602*10^-19;
X_b_Sb_0=0;
X_b_As_0=1;
X_s_Sb_0=0;
X_s_As_0=0;
temp=678;
phi_As=0.6*growth_rate;
phi_Sb=0.4*growth_rate;
v=10^13;
p1=v*exp(-E1/(kB*temp));
p2=v*exp(-E2/(kB*temp));
tspan = [0:0.01:1/growth_rate];
y0 = 0;
for c=0:ML+1
ML_counter=ML_counter+1;
Z=p1-p2;
L=-p1*X_s_As_0-2*p1*X_s_Sb_0-p1*X_b_Sb_0+p2*X_s_Sb_0-p2*X_b_As_0;
W=-p1*phi_As-2*p1*phi_Sb+p2*phi_Sb;
G=p1*X_s_Sb_0*phi_As+2*p1*X_s_Sb_0*phi_Sb+p1*X_s_As_0*phi_Sb+p1*X_b_Sb_0*phi_As+p1*X_b_Sb_0*phi_Sb;
T=p1*phi_Sb*phi_Sb+p1*phi_Sb*phi_As;
M=phi_Sb+p1*X_s_As_0*X_s_Sb_0+p1*X_s_Sb_0*X_s_Sb_0+p1*X_b_Sb_0*X_s_As_0+p1*X_b_Sb_0*X_s_Sb_0;
[t,y] = ode45(@(t,y) Z*y*y+L*y+W*y*t+G*t+T*t*t+M, tspan, y0);
XsIn=y(100/growth_rate+1);
XsGa=(phi_As+phi_Sb)*t(100/growth_rate+1)-XsIn;
XbGa=-XsGa+X_b_As_0+phi_As*t(100/growth_rate+1);
XbIn=-XsIn+X_b_Sb_0+phi_Sb*t(100/growth_rate+1);
X_b_Sb_0=XsIn;
X_b_As_0=XsGa;
X_b_Sb(1,ML_counter)=XbIn;
if (ML_counter>8)&&(ML_counter<18)
phi_As=1*growth_rate;
phi_Sb=0;    
end
discrete=discrete+(X_b_Sb(1,ML_counter)-X_b_Sb(1,ML_counter-1))*heaviside(time-ML_counter+3);
end
ml=[-1:1:ML+1];
plot(time,discrete,'m-','linewidth',1);
title('InAsSb/InAs: growth rate=1 : temp=420')
xlabel('Thickness(ML)')
ylabel('In concentration')
axis([-1 20 0 0.42])
hold on
growth_rate=1;
ML=20;
ML_counter=3;
discrete=0;
time=-1:0.001:20;
kB=1.38064852*10^-23;
E1=1.68*1.602*10^-19;
E2=1.75*1.602*10^-19;
X_b_Sb_0=0;
X_b_As_0=1;
X_s_Sb_0=0;
X_s_As_0=0;
temp=678;
phi_As=0.6*growth_rate;
phi_Sb=0.4*growth_rate;
v=10^13;
p1=v*exp(-E1/(kB*temp));
p2=v*exp(-E2/(kB*temp));
tspan = [0:0.01:1/growth_rate];
y0 = 0;
for c=0:ML+1
ML_counter=ML_counter+1;
Z=p1-p2;
L=-p1*X_s_As_0-2*p1*X_s_Sb_0-p1*X_b_Sb_0+p2*X_s_Sb_0-p2*X_b_As_0;
W=-p1*phi_As-2*p1*phi_Sb+p2*phi_Sb;
G=p1*X_s_Sb_0*phi_As+2*p1*X_s_Sb_0*phi_Sb+p1*X_s_As_0*phi_Sb+p1*X_b_Sb_0*phi_As+p1*X_b_Sb_0*phi_Sb;
T=p1*phi_Sb*phi_Sb+p1*phi_Sb*phi_As;
M=phi_Sb+p1*X_s_As_0*X_s_Sb_0+p1*X_s_Sb_0*X_s_Sb_0+p1*X_b_Sb_0*X_s_As_0+p1*X_b_Sb_0*X_s_Sb_0;
[t,y] = ode45(@(t,y) Z*y*y+L*y+W*y*t+G*t+T*t*t+M, tspan, y0);
XsIn=y(100/growth_rate+1);
XsGa=(phi_As+phi_Sb)*t(100/growth_rate+1)-XsIn;
XbGa=-XsGa+X_b_As_0+phi_As*t(100/growth_rate+1);
XbIn=-XsIn+X_b_Sb_0+phi_Sb*t(100/growth_rate+1);
X_b_Sb_0=XsIn;
X_b_As_0=XsGa;
X_b_Sb(1,ML_counter)=XbIn;
if (ML_counter>8)&&(ML_counter<18)
phi_As=1*growth_rate;
phi_Sb=0;    
end
discrete=discrete+(X_b_Sb(1,ML_counter)-X_b_Sb(1,ML_counter-1))*heaviside(time-ML_counter+3);
end
ml=[-1:1:ML+1];
plot(time,discrete,'b-','linewidth',1);
xlabel('Thickness(ML)')
ylabel('Sb concentration')
axis([-1 20 0 0.42])
hold on
R=0.6;
x0=0.4;
ML=60;
time=-20:0.001:40;
ML_counter=1;
X(1,1)=0;
for c=0:ML+1
    if ML_counter<7
    X(1,ML_counter+1)=(x0)*(1-R^(ML_counter));  
    end
    if 6<ML_counter
    X(1,ML_counter+1)=(x0)*(1-R^6)*R^(ML_counter-6);  
    end
    ML_counter=ML_counter+1;
end
plot(X,'k-','linewidth',1)
legend('kinetic:E1=1.6 ,E2=1.75','kinetic:E1=1.68 ,E2=1.75','Muraki:R=0.8')