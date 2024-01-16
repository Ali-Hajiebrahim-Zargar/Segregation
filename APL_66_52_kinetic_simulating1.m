clear all
clc;
% pattern:
% substrate-->GaAs
% 20 ML In0.2Ga0.8As
% growth rate=0.1 or 1 ML/s
% growth temperature = 400 or 500 c
% E1=1.8 ev
% E2=2 ev
%------  different temperature in growth rate=1   ----------
growth_rate=1;
ML=20;
ML_counter=2;                     
discrete=0;                       % to plot X_b_In discreetly
time=-1.01:0.001:25;
kB=1.38064852*10^-23;
E1=1.8*1.602*10^-19;
E2=2*1.602*10^-19;
X_b_In_0=0;                       %initial condition
X_b_Ga_0=1;                       %initial condition
X_s_In_0=0;                       %initial condition
X_s_Ga_0=0;                       %initial condition
temp=793.15;                      %growth temperature in kelvin
phi_Ga=0.8*growth_rate;           %effects of different growth rate
phi_In=0.2*growth_rate;           %effects of different growth rate
v=10^13;                          %Vibration Frequency for all semiconductors
p1=v*exp(-E1/(kB*temp));
p2=v*exp(-E2/(kB*temp));
tspan = [0:0.01:1/growth_rate];
y0 = 0;
for c=0:ML+1                      % for loop to calculate X_b_In 
ML_counter=ML_counter+1;
Z=p1-p2;
L=-p1*X_s_Ga_0-2*p1*X_s_In_0-p1*X_b_In_0+p2*X_s_In_0-p2*X_b_Ga_0;
W=-p1*phi_Ga-2*p1*phi_In+p2*phi_In;
G=p1*X_s_In_0*phi_Ga+2*p1*X_s_In_0*phi_In+p1*X_s_Ga_0*phi_In+p1*X_b_In_0*phi_Ga+p1*X_b_In_0*phi_In;
T=p1*phi_In*phi_In+p1*phi_In*phi_Ga;
M=phi_In+p1*X_s_Ga_0*X_s_In_0+p1*X_s_In_0*X_s_In_0+p1*X_b_In_0*X_s_Ga_0+p1*X_b_In_0*X_s_In_0;
[t,y] = ode45(@(t,y) Z*y*y+L*y+W*y*t+G*t+T*t*t+M, tspan, y0);
XsIn=y(100/growth_rate+1);                           %XsIn in t=1
XsGa=(phi_Ga+phi_In)*t(100/growth_rate+1)-XsIn;      %use kinetic's equations to calculate next ML's initial condition
XbGa=-XsGa+X_b_Ga_0+phi_Ga*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
XbIn=-XsIn+X_b_In_0+phi_In*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
X_b_In_0=XsIn;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_Ga_0=XsGa;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_In(1,ML_counter)=XbIn;                           %save results 
if ML_counter>ML                                     % change In and Ga phi if 20 ML finished
phi_Ga=0;
phi_In=0;    
end
discrete=discrete+(X_b_In(1,ML_counter)-X_b_In(1,ML_counter-1))*heaviside(time-ML_counter+3);   % use stepp function to plot discreetly
end
ml=[-1:1:ML+1];                                       %print X_b_In for each monolayer
figure
plot(time,discrete,'m-','linewidth',1)
title('In0.2Ga0.8As/GaAs growth rate=1: temp=500 & temp=400')
xlabel('Thickness(ML)')
ylabel('In concentration')
axis([-1 21 0 1])
hold on
R=0.75;
x0=0.2;
ML=60;
time=-20:0.001:40;
ML_counter=1;
X(1,1)=0;
for c=0:ML+1
    X(1,ML_counter+1)=(x0)*(1-R^(ML_counter));  
    ML_counter=ML_counter+1;
end
ml=[-1:1:ML+1];
plot(X,'b-','linewidth',1)
title('In0.36 Ga0.64 N0.04 As0.96')
xlabel('Thickness(ML)')
ylabel('In concentration')
%end of plotting for temp=500 c and growth rate=1
growth_rate=1;
ML=20;
ML_counter=2;                     
discrete=0;                       % to plot X_b_In discreetly
time=-1.01:0.001:25;
kB=1.38064852*10^-23;
E1=1.8*1.602*10^-19;
E2=2*1.602*10^-19;
X_b_In_0=0;                       %initial condition
X_b_Ga_0=1;                       %initial condition
X_s_In_0=0;                       %initial condition
X_s_Ga_0=0;                       %initial condition
temp=673.15;                      %growth temperature in kelvin
phi_Ga=0.8*growth_rate;           %effects of different growth rate
phi_In=0.2*growth_rate;           %effects of different growth rate
v=10^13;                          %Vibration Frequency for all semiconductors
p1=v*exp(-E1/(kB*temp));
p2=v*exp(-E2/(kB*temp));
tspan = [0:0.01:1/growth_rate];
y0 = 0;
for c=0:ML+1                      % for loop to calculate X_b_In 
ML_counter=ML_counter+1;
Z=p1-p2;
L=-p1*X_s_Ga_0-2*p1*X_s_In_0-p1*X_b_In_0+p2*X_s_In_0-p2*X_b_Ga_0;
W=-p1*phi_Ga-2*p1*phi_In+p2*phi_In;
G=p1*X_s_In_0*phi_Ga+2*p1*X_s_In_0*phi_In+p1*X_s_Ga_0*phi_In+p1*X_b_In_0*phi_Ga+p1*X_b_In_0*phi_In;
T=p1*phi_In*phi_In+p1*phi_In*phi_Ga;
M=phi_In+p1*X_s_Ga_0*X_s_In_0+p1*X_s_In_0*X_s_In_0+p1*X_b_In_0*X_s_Ga_0+p1*X_b_In_0*X_s_In_0;
[t,y] = ode45(@(t,y) Z*y*y+L*y+W*y*t+G*t+T*t*t+M, tspan, y0);
XsIn=y(100/growth_rate+1);                           %XsIn in t=1
XsGa=(phi_Ga+phi_In)*t(100/growth_rate+1)-XsIn;      %use kinetic's equations to calculate next ML's initial condition
XbGa=-XsGa+X_b_Ga_0+phi_Ga*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
XbIn=-XsIn+X_b_In_0+phi_In*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
X_b_In_0=XsIn;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_Ga_0=XsGa;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_In(1,ML_counter)=XbIn;                           %save results 
if ML_counter>ML                                     % change In and Ga phi if 20 ML finished
phi_Ga=0;
phi_In=0;    
end
discrete=discrete+(X_b_In(1,ML_counter)-X_b_In(1,ML_counter-1))*heaviside(time-ML_counter+3);   % use stepp function to plot discreetly
end
ml=[-1:1:ML+1];
plot(time,discrete,'b-','linewidth',1)
legend('temp=500','temp=400')
hold off
%end of plotting for temp=400 c and growth rate=1
%--------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------
%----------------------  different growth rate in temp=500   --------------------------------

growth_rate=1;
ML=20;
ML_counter=1;                     
discrete=0;                       % to plot X_b_In discreetly
time=-1.01:0.001:25;
kB=1.38064852*10^-23;
E1=1.8*1.602*10^-19;
E2=2*1.602*10^-19;
X_b_In_0=0;                       %initial condition
X_b_Ga_0=1;                       %initial condition
X_s_In_0=0;                       %initial condition
X_s_Ga_0=0;                       %initial condition
temp=773.15;                      %growth temperature in kelvin
phi_Ga=0.8*growth_rate;           %effects of different growth rate
phi_In=0.2*growth_rate;           %effects of different growth rate
v=10^13;                          %Vibration Frequency for all semiconductors
p1=v*exp(-E1/(kB*temp));
p2=v*exp(-E2/(kB*temp));
tspan = [0:0.01:1/growth_rate];
y0 = 0;
for c=0:ML+1                      % for loop to calculate X_b_In 
ML_counter=ML_counter+1;
Z=p1-p2;
L=-p1*X_s_Ga_0-2*p1*X_s_In_0-p1*X_b_In_0+p2*X_s_In_0-p2*X_b_Ga_0;
W=-p1*phi_Ga-2*p1*phi_In+p2*phi_In;
G=p1*X_s_In_0*phi_Ga+2*p1*X_s_In_0*phi_In+p1*X_s_Ga_0*phi_In+p1*X_b_In_0*phi_Ga+p1*X_b_In_0*phi_In;
T=p1*phi_In*phi_In+p1*phi_In*phi_Ga;
M=phi_In+p1*X_s_Ga_0*X_s_In_0+p1*X_s_In_0*X_s_In_0+p1*X_b_In_0*X_s_Ga_0+p1*X_b_In_0*X_s_In_0;
[t,y] = ode45(@(t,y) Z*y*y+L*y+W*y*t+G*t+T*t*t+M, tspan, y0);
XsIn=y(100/growth_rate+1);                           %XsIn in t=1
XsGa=(phi_Ga+phi_In)*t(100/growth_rate+1)-XsIn;      %use kinetic's equations to calculate next ML's initial condition
XbGa=-XsGa+X_b_Ga_0+phi_Ga*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
XbIn=-XsIn+X_b_In_0+phi_In*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
X_b_In_0=XsIn;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_Ga_0=XsGa;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_In(1,ML_counter)=XbIn;                           %save results 
if ML_counter>ML                                     % change In and Ga phi if 20 ML finished
phi_Ga=0;
phi_In=0;    
end
discrete=discrete+(X_b_In(1,ML_counter)-X_b_In(1,ML_counter-1))*heaviside(time-ML_counter+3);   % use stepp function to plot discreetly
end
ml=[-1:1:ML+1];
[ml;X_b_In]
figure
plot(time,discrete,'m-','linewidth',1)
title('In0.2Ga0.8As/GaAs temp=500: growth rate=1 & growth rate=0.1')
xlabel('Thickness(ML)')
ylabel('In concentration')
axis([-1 21 0 1])
hold on
%end of plotting for temp=500 c and growth rate=1
growth_rate=0.1;
ML=20;
ML_counter=1;                     
discrete=0;                       % to plot X_b_In discreetly
time=-1.01:0.001:25;
kB=1.38064852*10^-23;
E1=1.8*1.602*10^-19;
E2=2*1.602*10^-19;
X_b_In_0=0;                       %initial condition
X_b_Ga_0=1;                       %initial condition
X_s_In_0=0;                       %initial condition
X_s_Ga_0=0;                       %initial condition
temp=773.15;                      %growth temperature in kelvin
phi_Ga=0.8*growth_rate;           %effects of different growth rate
phi_In=0.2*growth_rate;           %effects of different growth rate
v=10^13;                          %Vibration Frequency for all semiconductors
p1=v*exp(-E1/(kB*temp));
p2=v*exp(-E2/(kB*temp));
tspan = [0:0.01:1/growth_rate];
y0 = 0;
for c=0:ML+1                      % for loop to calculate X_b_In 
ML_counter=ML_counter+1;
Z=p1-p2;
L=-p1*X_s_Ga_0-2*p1*X_s_In_0-p1*X_b_In_0+p2*X_s_In_0-p2*X_b_Ga_0;
W=-p1*phi_Ga-2*p1*phi_In+p2*phi_In;
G=p1*X_s_In_0*phi_Ga+2*p1*X_s_In_0*phi_In+p1*X_s_Ga_0*phi_In+p1*X_b_In_0*phi_Ga+p1*X_b_In_0*phi_In;
T=p1*phi_In*phi_In+p1*phi_In*phi_Ga;
M=phi_In+p1*X_s_Ga_0*X_s_In_0+p1*X_s_In_0*X_s_In_0+p1*X_b_In_0*X_s_Ga_0+p1*X_b_In_0*X_s_In_0;
[t,y] = ode45(@(t,y) Z*y*y+L*y+W*y*t+G*t+T*t*t+M, tspan, y0);
XsIn=y(100/growth_rate+1);                           %XsIn in t=1
XsGa=(phi_Ga+phi_In)*t(100/growth_rate+1)-XsIn;      %use kinetic's equations to calculate next ML's initial condition
XbGa=-XsGa+X_b_Ga_0+phi_Ga*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
XbIn=-XsIn+X_b_In_0+phi_In*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
X_b_In_0=XsIn;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_Ga_0=XsGa;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_In(1,ML_counter)=XbIn;                           %save results 
if ML_counter>ML                                     % change In and Ga phi if 20 ML finished
phi_Ga=0;
phi_In=0;    
end
discrete=discrete+(X_b_In(1,ML_counter)-X_b_In(1,ML_counter-1))*heaviside(time-ML_counter+3);   % use stepp function to plot discreetly
end
ml=[-1:1:ML+1];
[ml;X_b_In]
plot(time,discrete,'b-','linewidth',1)
legend('growth rate=1','growth rate=0.1')
hold off
%end of plotting for temp=500 c and growth rate=0.1
%--------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------
%----------------------  different growth rate in temp=400   --------------------------------
growth_rate=1;
ML=20;
ML_counter=1;                     
discrete=0;                       % to plot X_b_In discreetly
time=-1.01:0.001:25;
kB=1.38064852*10^-23;
E1=1.8*1.602*10^-19;
E2=2*1.602*10^-19;
X_b_In_0=0;                       %initial condition
X_b_Ga_0=1;                       %initial condition
X_s_In_0=0;                       %initial condition
X_s_Ga_0=0;                       %initial condition
temp=673.15;                      %growth temperature in kelvin
phi_Ga=0.8*growth_rate;           %effects of different growth rate
phi_In=0.2*growth_rate;           %effects of different growth rate
v=10^13;                          %Vibration Frequency for all semiconductors
p1=v*exp(-E1/(kB*temp));
p2=v*exp(-E2/(kB*temp));
tspan = [0:0.01:1/growth_rate];
y0 = 0;
for c=0:ML+1                      % for loop to calculate X_b_In 
ML_counter=ML_counter+1;
Z=p1-p2;
L=-p1*X_s_Ga_0-2*p1*X_s_In_0-p1*X_b_In_0+p2*X_s_In_0-p2*X_b_Ga_0;
W=-p1*phi_Ga-2*p1*phi_In+p2*phi_In;
G=p1*X_s_In_0*phi_Ga+2*p1*X_s_In_0*phi_In+p1*X_s_Ga_0*phi_In+p1*X_b_In_0*phi_Ga+p1*X_b_In_0*phi_In;
T=p1*phi_In*phi_In+p1*phi_In*phi_Ga;
M=phi_In+p1*X_s_Ga_0*X_s_In_0+p1*X_s_In_0*X_s_In_0+p1*X_b_In_0*X_s_Ga_0+p1*X_b_In_0*X_s_In_0;
[t,y] = ode45(@(t,y) Z*y*y+L*y+W*y*t+G*t+T*t*t+M, tspan, y0);
XsIn=y(100/growth_rate+1);                           %XsIn in t=1
XsGa=(phi_Ga+phi_In)*t(100/growth_rate+1)-XsIn;      %use kinetic's equations to calculate next ML's initial condition
XbGa=-XsGa+X_b_Ga_0+phi_Ga*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
XbIn=-XsIn+X_b_In_0+phi_In*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
X_b_In_0=XsIn;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_Ga_0=XsGa;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_In(1,ML_counter)=XbIn;                           %save results 
if ML_counter>ML                                     % change In and Ga phi if 20 ML finished
phi_Ga=0;
phi_In=0;    
end
discrete=discrete+(X_b_In(1,ML_counter)-X_b_In(1,ML_counter-1))*heaviside(time-ML_counter+3);   % use stepp function to plot discreetly
end
ml=[-1:1:ML+1];
[ml;X_b_In]
figure
plot(time,discrete,'m-','linewidth',1)
title('In0.2Ga0.8As/GaAs temp=400: growth rate=1 & growth rate=0.1')
xlabel('Thickness(ML)')
ylabel('In concentration')
axis([-1 21 0 1])
hold on
%end of plotting for temp=400 c and growth rate=1
growth_rate=0.1;
ML=20;
ML_counter=1;                     
discrete=0;                       % to plot X_b_In discreetly
time=-1.01:0.001:25;
kB=1.38064852*10^-23;
E1=1.8*1.602*10^-19;
E2=2*1.602*10^-19;
X_b_In_0=0;                       %initial condition
X_b_Ga_0=1;                       %initial condition
X_s_In_0=0;                       %initial condition
X_s_Ga_0=0;                       %initial condition
temp=673.15;                      %growth temperature in kelvin
phi_Ga=0.8*growth_rate;           %effects of different growth rate
phi_In=0.2*growth_rate;           %effects of different growth rate
v=10^13;                          %Vibration Frequency for all semiconductors
p1=v*exp(-E1/(kB*temp));
p2=v*exp(-E2/(kB*temp));
tspan = [0:0.01:1/growth_rate];
y0 = 0;
for c=0:ML+1                      % for loop to calculate X_b_In 
ML_counter=ML_counter+1;
Z=p1-p2;
L=-p1*X_s_Ga_0-2*p1*X_s_In_0-p1*X_b_In_0+p2*X_s_In_0-p2*X_b_Ga_0;
W=-p1*phi_Ga-2*p1*phi_In+p2*phi_In;
G=p1*X_s_In_0*phi_Ga+2*p1*X_s_In_0*phi_In+p1*X_s_Ga_0*phi_In+p1*X_b_In_0*phi_Ga+p1*X_b_In_0*phi_In;
T=p1*phi_In*phi_In+p1*phi_In*phi_Ga;
M=phi_In+p1*X_s_Ga_0*X_s_In_0+p1*X_s_In_0*X_s_In_0+p1*X_b_In_0*X_s_Ga_0+p1*X_b_In_0*X_s_In_0;
[t,y] = ode45(@(t,y) Z*y*y+L*y+W*y*t+G*t+T*t*t+M, tspan, y0);
XsIn=y(100/growth_rate+1);                           %XsIn in t=1
XsGa=(phi_Ga+phi_In)*t(100/growth_rate+1)-XsIn;      %use kinetic's equations to calculate next ML's initial condition
XbGa=-XsGa+X_b_Ga_0+phi_Ga*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
XbIn=-XsIn+X_b_In_0+phi_In*t(100/growth_rate+1);     %use kinetic's equations to calculate next ML's initial condition
X_b_In_0=XsIn;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_Ga_0=XsGa;                                       %use kinetic's equations to calculate next ML's initial condition
X_b_In(1,ML_counter)=XbIn;                           %save results 
if ML_counter>ML                                     % change In and Ga phi if 20 ML finished
phi_Ga=0;
phi_In=0;    
end
discrete=discrete+(X_b_In(1,ML_counter)-X_b_In(1,ML_counter-1))*heaviside(time-ML_counter+3);   % use stepp function to plot discreetly
end
ml=[-1:1:ML+1];
[ml;X_b_In]
plot(time,discrete,'b-','linewidth',1)
legend('growth rate=1','growth rate=0.1')
hold off
%end of plotting for temp=400 c and growth rate=0.1
%--------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------
