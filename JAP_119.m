clc
clear all
R=0.81;
x0=0.4;
ML=70;
ML_counter=0;
X(1,1)=0;
for c=0:ML
    if ML_counter<7
    X(1,ML_counter+1)=(x0)*(1-R^(ML_counter));  
    end
    if (6<ML_counter)&&(ML_counter<18)
    X(1,ML_counter+1)=(x0)*(1-R^6)*R^(ML_counter-6);  
    end
    if (17<ML_counter)&&(ML_counter<24)
    X(1,ML_counter+1)=(x0)*(1-R^(ML_counter-17)+((x0)*(1-R^6)*R^11)*R^(ML_counter-17));  
    end
    if (23<ML_counter)&&(ML_counter<35)
    X(1,ML_counter+1)=(x0)*(1-R^(23-17)+((x0)*(1-R^6)*R^11)*R^(23-17))*R^(ML_counter-23);  
    end
    if (34<ML_counter)&&(ML_counter<41)
    X(1,ML_counter+1)=(x0)*(1-R^(ML_counter-34));  
    end
    if (40<ML_counter)&&(ML_counter<52)
    X(1,ML_counter+1)=(x0)*(1-R^40)*R^(ML_counter-40);  
    end
     if (51<ML_counter)&&(ML_counter<58)
    X(1,ML_counter+1)=(x0)*(1-R^(ML_counter-51));  
    end
    if (57<ML_counter)&&(ML_counter<69)
    X(1,ML_counter+1)=(x0)*(1-R^57)*R^(ML_counter-57);  
    end
    ML_counter=ML_counter+1;
end
X
figure
plot(X)
hold on
