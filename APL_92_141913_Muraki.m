clc
clear all
R=0.68;
x0=0.36;
ML=60;
time=-20:0.001:40;
ML_counter=1;
X(1,1)=0;
for c=0:ML+1
    if ML_counter<17
    X(1,ML_counter+1)=(x0)*(1-R^(ML_counter));  
    end
    if 16<ML_counter
    X(1,ML_counter+1)=(x0)*(1-R^16)*R^(ML_counter-16);  
    end
    ML_counter=ML_counter+1;
end
ml=[-1:1:ML+1];
[ml;X]
figure
plot(X,'m-','linewidth',1)
title('In0.36 Ga0.64 N0.04 As0.96')
xlabel('Thickness(ML)')
ylabel('In concentration')
axis([-20 45 -0.05 0.4])
