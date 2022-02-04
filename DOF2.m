%% 2DOF micro Model PZT-5A

% Secondary beam
% Set the dimentions and determine the mass then the stiffness
% Using PMMA and PZT
clear all
close all
clc
E = 3200000000;
Y = 5.4e10;
a = 9.81;
t2 = 0.00001;
w22 = 1.3*(2*pi); % frequency of the secondary beam, equal to the main frequency of the system
s2 = 8150000;
l2 = 0.022;
b2 = 0.001
lp = l2/2;
bp = b2 - 1e-4;
tp = 10e-6;
In1 = (b2*t2^3)/12;
In2 = (bp*tp^3)/12;
Db = 1190;
Dp = 1780;
mb = Db*l2*b2*t2;
mp = Dp*lp*bp*tp;
kb = (3*E*In1)/(l2^3);
kp = (3*Y*In2)/(lp^3);
k2 = 1/((1/kb)+(1/kp))
%k2 = kb
m2 = 6e-6;% k2/(w22^2)%-0.25*mb
%m2 = (s2*b2*(t2^2))/(6*a*l2)
%k2 = (m2*(w22^2))
Vb2 = l2*b2*t2;
Vp = lp*bp*tp;

s2 = (6*m2*a*l2)/(b2*(t2^2))

% Mass_2 dimentions
D = 8960;
v2 = m2/D;
W2 = b2;
H2 = 0.00004;
L2 = v2/(W2*H2);
%%

% Primary beam
% Determine the mass and stiffness
syms m1 k1
w1 = (1.3 - (1.3*0.2))*(2*pi);
w2 = (1.3 + (1.3*0.2))*(2*pi);
t1 = t2;
Eq1 = m1*(m2*(w1^2)*(w2^2)) - k1*k2 == 0;
Eq2 = m1*(((w1^2)+(w2^2))*m2 - k2) - k1*m2 - m2*k2 == 0;
Sol = solve(Eq1, Eq2, m1, k1);
k1 = vpa(Sol.k1)
m1 = vpa(Sol.m1)
% Set the thickness and length then determine the width
l1 = 0.026;
b1 = vpa((4*(l1^3)*k1)/(E*(t1^3)))
s1 = vpa((6*m1*a*l1)/(b1*(t1^2)))
bb1 = vpa((4*(l1^3)*(k1/2))/(E*(t1^3)))
Vb1 = l1*t2*b1;

% Dimension of mass 1
v1 = m1/D;
W1 = b1;
H1 = 0.00002;
L1 = v1/(W1*H1);
Vt = Vb1+Vb2+v1+v2+(3*Vp);

%%
%Mecanical Power Calculation
W = [0:0.1:2];
Ww = (2*pi)*W;
n = length(Ww);
A = a*ones(1,n);
for i=1:n
    X(i) = (-1*A(i))/((Ww(i))^2);
end

for i=1:n
    Xx(i) = (k1+k2-m1*(Ww(i))^2)*(k2-m2*(Ww(i))^2)-(k2)^2;
    X1(i) = ((k2-m2*(Ww(i))^2)*k1*X(i))/Xx(i);
end

for i=1:n
    Xx(i) = (k1+k2-m1*(Ww(i))^2)*(k2-m2*(Ww(i))^2)-(k2)^2;
    X2(i) = (k1*k2*X(i))/Xx(i);
end

Fig1 = figure('Name','Amplitude')
plot(W,X,W,X1,'r--',W,X2,'k:','LineWidth',2)
xlabel('Frequency(Hz)')
ylabel('Amplitude')
legend('Input','Primary beam','Secondary beam','Location','best')
axis([0.8 1.8 0 5])

Dis1 = X-X1;
Dis2 = X2-X1;
Fig = figure('Name','Displacement')
plot(W,Dis1,'r--',W,Dis2,'k:','LineWidth',2)
xlabel('Frequency(Hz)')
ylabel('Displacement (m)')
legend('X-X1','X2-X1')
axis([0.8 1.8 -4 3])

P1 = ones(n,1);
P2 = ones(n,1);
c = 0.1;
for i=1:n
    P1(i) = c*((Ww(i))^2)*((X(i)-X1(i))^2);
    P2(i) = c*((Ww(i))^2)*((X2(i)-X1(i))^2);
end

Fig2 = figure ('Name','Mechanical Power')
plot(W,P1,W,P2,'r--','LineWidth',2)
xlabel('Frequency(Hz)')
ylabel('Power(W)')
legend('Primary beam','Secondary beam')

%%
%Voltage calculation
d1 = X-X1;
d31 = -190e-12;
K = 0.34;
Kt = 1800;
e = 8.854e-12; %piezoelectric permitivity
e33 = Kt*e; %transvers dielectric constant
Ca = (e33*lp*bp)/tp;% (((lp/0.0254)*(bp/0.0254)*Kt)/(4.45*(tp/0.0254)))*1e-9;
t = t2;
for i=1:n
    F1(i)= (E*b1*t^3*d1(i))/(4*l1^3);
end
for i=1:n
    sp(i)= (6*F1(i)*l1)/(b1*t^2);
end

%Y = 4e+9;
Voc1 = ones(n,1);
for i=1:n
    a1 = -3/(1-K^2);
    a2 = (d31*Y*tp)/(e33*l1);
    a3 = ((l1+L1)*t)/(4*(l1^2)+9*l1*L1+6*L1^2);
    Voc1(i) = 2*a1*a2*a3*d1(i);
end


d2 = X2-X1;
for i=1:n
    F2(i)= (E*b2*t^3*d2(i))/(4*l2^3);
end
for i=1:n
    sp(i)= (6*F1(i)*l2)/(b2*t^2);
end

Voc2 = ones(n,1);
for i=1:n
    a1 = -3/(1-K^2);
    a2 = (d31*Y*tp)/(e33*l2);
    a3 = ((l2+L2)*t)/(4*l2^2+9*l2*L2+6*L2^2);
    Voc2(i) = a1*a2*a3*d2(i);
end

VocT = Voc2-Voc1;
Voc = sum(VocT(10:18))/9
PMc = P1+P2;
Fig3 = figure('Name','Open circuit voltage')
plot(W,Voc1,W,Voc2,'r--','LineWidth',2)
xlabel('Frequency(Hz)')
ylabel('Voltage(V)')
legend('Primary beam','Secondary beam')

%%
%Electrical Power Calculation
Re = 1/(Ca*w22*2*pi);
for i=1:n
    Z(i) = 1/(Ca*Ww(i));
end
RL1 = 900;
Rt1 = RL1*ones(1,n)+Z;

V1 = ones(n,1);
Po1 = ones(n,1);
for i=1:n
    In1(i) = Voc1(i)/Rt1(i);
    V1(i) = In1(i)*RL1;
    Po1(i) = ((Voc1(i)^2)*RL1)/(2*(Rt1(i))^2);
end

V2 = ones(n,1);
Po2 = ones(n,1);
for i=1:n
    I2(i) = Voc2(i)/Rt1(i);
    V2(i) = I2(i)*RL1;
   Po2(i) = ((Voc2(i)^2)*RL1)/(2*(Rt1(i))^2);
end

VT1 = V1 + V2;
Vmoy = sum(VT1(10:18))/8
Fig4 = figure('Name','Output voltage')
plot(W,VT1,W,V1,'r--',W,V2,'k:','LineWidth',2)
xlabel('Frequency(Hz)')
ylabel('Voltage(V)')
legend('Total','Primary beam','Secondary beam')

PT1 = Po1 + Po2;
moy = PT1(10:18);
Pmoy = sum(moy)/9
Fig5 = figure('Name','Electrical Power')
plot(W,PT1,W,Po1,'r--',W,Po2,'k:','LineWidth',2)
xlabel('Frequency(Hz)')
ylabel('Power(W)')
legend('Total','Primary beam','Secondary beam')

%Maximum volyage and power
RL2 = 2000*ones(1,n);
Rt2 = RL2+Z;

V12 = ones(n,1);
Po12 = ones(n,1);
for i=1:n
    In1(i) = Voc1(i)/Rt2(i);
    V12(i) = In1(i)*RL2(i);
    Po12(i) = ((V12(i)^2)*RL2(i))/(2*(Rt2(i))^2);
end

V22 = ones(n,1);
Po22 = ones(n,1);
for i=1:n
    I2(i) = Voc2(i)/Rt2(i);
    V22(i) = I2(i)*RL2(i);
   Po22(i) = ((V22(i)^2)*RL2(i))/(2*(Rt2(i))^2);
end

VT2 = V12 + V22;
%Fig6 = figure('Name','Maximum Output voltage')
%plot(W,VT2,W,V12,'r--',W,V22,'k:','LineWidth',2)
%xlabel('Frequency(Hz)')
%ylabel('Voltage(V)')
%legend('Total','Primary beam','Secondary beam')

PT2 = Po12 + Po22;
%Fig7 = figure('Name','Maximum Electrical Power')
%plot(W,PT2,W,Po12,'r--',W,Po22,'k:','LineWidth',2)
%xlabel('Frequency(Hz)')
%ylabel('Power(W)')
%legend('Total','Primary beam','Secondary beam')


%Power density

PD1 = ones(n,1);
PD2 = ones(n,1);
for i=1:n
    PD1(i) = (1e6)*((PT1(i))/Vt);
    PD2(i) = 1e6*((PT2(i))/Vt);
end
F8 = figure('Name','Power density')
plot(W,PD1)
xlabel('Frequency(Hz)')
ylabel('Power density(uW/cm3)')

%F9 = figure('Name','Maximum Power density')
%plot(W,PD2)
%xlabel('Frequency(Hz)')
%ylabel('Power density(uW/cm3)')

%Energy Calculation

E11 = ones(n,1);
E2 = ones(n,1);
E = ones(n,1);
for i=1:n
    E11(i) = (Po1(i))/(Ww(i));
    E2(i) = (Po2(i))/(Ww(i));
    Ep1(i) = (PT1(i))/(Ww(i));
    Ep2(i) = (PT2(i))/(Ww(i));
end
ET1 = sum(Ep1);
ET2 = sum(Ep2);
%Fig10 = figure('Name','Energy')
%plot(W,Ep1,W,Ep2,'r',W,E11,'--',W,E2,'k:','LineWidth',2)
%xlabel('Frequency(Hz)')
%ylabel('Energy(J)')
%legend('Total','Maximum','Primary beam','Secondary beam')

ED1 = ones(n,1);
ED2 = ones(n,1);
for i=1:n
    ED1(i) = Ep1(i)/Vt;
    ED2(i) = Ep2(i)/Vt;
end

%Fig11 = figure('Name','Energy density')
%plot(W,ED1,W,ED2,'--')
%xlabel('Frequency(Hz)')
%ylabel('Energy density (J/m3)')
%legend('Normal','Maximu')

xlswrite('2DOF_PZT_Voltages',[V1(:),V2(:),VT1(:),Voc1(:),Voc2(:),VocT(:)])
xlswrite('2DOF_PZT_Power_Output',[Po1(:),Po2(:),PT1(:),PD1(:),ED1(:)])
xlswrite('2DOF_PZT_Power_Mecha',[P1(:),P2(:)])