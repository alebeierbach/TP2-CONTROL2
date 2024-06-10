clear all; clc; close all;

filename='C:\Users\Ale\Desktop\CONTROL 2\TP1 CONTROL2\Curvas_Medidas_Motor_2024.xls';
sheet=1;
total=xlsread(filename,sheet,'A1:E32370');
tdata=total(:,1);
wdata=total(:,2);
idata=total(:,3);
vdata=total(:,4);
tldata=total(:,5);
%plot(tdata,tldata)

s=tf('s');
Ra=55.6;
La=366e-6;
Km=6.53e-3;
Ki=6.49e-3;
Bm=0;
J=5e-9;

A=[-Ra/La -Km/La 0;
    Ki/J -Bm/J 0;
    0 1 0];
B=[1/La 0 ;
    0 -1/J;
    0 0];
C= [0 0 1];
D= [0 0];

sys=ss(A,B,C,D);
FdT=tf(sys)

%Defino las matrices ampliadas:
Aa=[A zeros(3,1);
    -C 0];
Ba=[B(:,1);0]; %Saco al Tl como entrada, ahora sera una perturbacion.
Ca=[C 0];

%Defino Q y R para LQR:
Q=diag([10 10 10 100000000]);
R=1;
K=lqr(Aa,Ba,Q,R)    %Calculo controlador K con LQR e integrador

%%
%Ahora calculare el tiempo de integracion de euler:
polos=eig(Aa-Ba*K);
lambda=max(polos);  %Hallo el valor maximo de polo
%tiempo de integracion: Se simulara con un tiempo menor al mismo
tr=log(0.95)/lambda;

h=4e-07;
t=0:h:0.6-h;   %Defino mi vector de tiempos para la simulacion

%Defino mi entrada referencia del sistema..
reference= pi/2*square(2*pi*1/0.6*t); %T=1s
plot(t,reference,'LineWidth',1.5)
xlabel('Tiempo [seg]')
ylabel('Ángulo [rad]')
title('Señal de referencia \theta_i')
grid
%%
%Simulacion del controlador

ia(1)=0;
w(1)=0;
tita(1)=0;

x=[ia(1) w(1) tita(1)]; %vector de estados auxiliar
statevector=[ia(1) w(1) tita(1)]'; %vector de estados xt
dX=[0 0 0];
zeta(1) = 0;
integ(1) = zeta(1);
k=1;

for i=1:(0.6/h)
    zetaP = reference(i)-C*statevector; %Integrador
    zeta(i) = integ+zetaP*h; %Mi integrador se va actualizando adecuando la salida de control
    u(i)=-K(1:3)*statevector-K(4)*zeta(i);    %accion de control
    ia(i)=x(1);
    w(i)=x(2);
    tita(i)=x(3);
    x1p = -Ra*x(1)/La-Km*x(2)/La+u(i)/La;    %se aplica accion de control a mi sist.
    x2p = Ki*x(1)/J-Bm*x(2)/J-tldata(k)/J;        %resto el torque
    x3p = x(2);
    dX=[x1p x2p x3p];
    x=x+h*dX;   %Voy actualizando mi vector de estados
    statevector=[ia(i) w(i) tita(i)]';
    integ = zeta(i);
    k=k+1;
    if(k==32369)    %Me aseguro que el torque este en bucle..
        k=1;
    end
end

plot(t,tita);
xlabel('Tiempo [seg]')
title('Angulo')
grid
hold on;
plot(t,u,'LineWidth',1.5)
xlabel('Tiempo [seg]')
title('Acción de control')
grid
hold off;

%% Sistema con observador:

Ra=55.6;
Laa=366e-6;
Km=6.53e-3;
Ki=6.49e-3;
Bm=0;
J=5e-9;
A = [-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B = [1/Laa 0 0]';
C = [0 0 1];
D = [0 0];
Q = diag([0.1 0.1 0.1 100000]);
R = 100;
Aamp = [A zeros(3,1); -C 0];
Bamp = [B(:,1); 0];
Camp = [C 0];
K = lqr(Aamp,Bamp,Q,R);
reference = (pi/2)*square(2*pi*(1/50)*t);
torque = ((1.15e-3)/2)*square(2*pi*(1/50)*t)+((1.15e-3)/2);
ia(1) = 0;
theta(1) = 0;
omega(1) = 0;
stateVector = [ia(1) omega(1) theta(1)]';
xop = [0 0 0]';
x = [ia(1) omega(1) theta(1)]';
zeta(1) = 0;
integ(1) = zeta(1);
%% 
% Ahora, para el observador:

Ao = A'                      %Todos las matrices tienen que ser convertidas
Bo = C'
Co = B' 
Qo = diag([.1 .1 .1])
Ro = 10;
Ko = lqr(Ao,Bo,Qo,Ro)   %Como para el observador no implementamos integrador, quedara de orden 3
obsStateVector = [ia(1) omega(1) theta(1)]';    %seria mi xhat, vector de estados construido en base a observar la salida
xObs = [0 0 0]';    %Aca se guardaran los valores observados


for i=1:(0.6/h)
    zetaP = reference(i)-Camp(1:3)*stateVector-Camp(4)*integ;
    zeta(i) = integ+zetaP*h;
    u(i) = -K(1:3)*obsStateVector-K(4)*zeta(i);     %utilizo el vector de estados observado
    ia(i) = x(1);
    omega(i) = x(2);
    theta(i) = x(3);
    x1P = -Ra*x(1)/Laa-Km*x(2)/Laa+u(i)/Laa;
    x2P = Ki*x(1)/J-Bm*x(2)/J-torque(i)/J;
    x3P = x(2);
    xP = [x1P x2P x3P]';
    x = x+xP*h;
    iaO(i)= xObs(1);
    omegaO(i)= xObs(2);
    thetaO(i)= xObs(3);
    yO(i) = C*obsStateVector;
    y(i) = Camp(1:3)*stateVector+Camp(4)*integ;
    xTP = A*xObs+B*u(i)+Ko*(y(:,i)-yO(:,i));
    xObs = xObs+xTP*h;
    stateVector = [ia(i) omega(i) theta(i)]';
    integ = zeta(i);
    obsStateVector =[iaO(i) omegaO(i) thetaO(i)]';
    k=k+1;
    if(k==32369)    %Me aseguro que el torque este en bucle..
        k=1;
    end
end

plot(t,iaO)     %Grafico de la corriente observada






