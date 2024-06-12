% Etapa II: Definición del Modelo de Estado para el Sistema Continuo

%-------------------------------------------------------------------------%
%------------------------- Parámetros del Avión --------------------------%
%-------------------------------------------------------------------------%
param_a = 0.07;
param_b = 5;
velocidad_vuelo = 150;    % [m/s]
frecuencia_natural = 9;   % [rad/s]

%-------------------------------------------------------------------------%
%----------------------- Modelado en Espacio de Estados ------------------%
%-------------------------------------------------------------------------%
% Condiciones iniciales de las variables de estado
angulo_alpha(1) = 0;
angulo_phi(1)   = 0;
vel_phi(1)      = 0;
altura_h(1)     = 500;
estado_inicial = [angulo_alpha(1) angulo_phi(1) vel_phi(1) altura_h(1)]';

% Matriz del Sistema
matriz_Ac = [-param_a    param_a   0 0;    
              0           0        1 0;    
              frecuencia_natural^2 -frecuencia_natural^2  0 0;    
              velocidad_vuelo     0        0 0];

% Matriz de Entrada
matriz_Bc = [0 0 param_b*frecuencia_natural^2 0]';

% Matriz de Salida
matriz_Cc = [0 0 0 1;           
              0 1 0 0];         

% Matriz de Acoplamiento Directo
matriz_Dc = 0;

%-------------------------------------------------------------------------%
%----------------- Representación del Sistema en Matlab ------------------%
%-------------------------------------------------------------------------%
sistema_continuo  = ss(matriz_Ac, matriz_Bc, matriz_Cc, matriz_Dc);     

%%
% Etapa III: Análisis de Controlabilidad

%-------------------------------------------------------------------------%
%------------- Verificación de Controlabilidad del Sistema ---------------%
%-------------------------------------------------------------------------%

% Comprobación de Controlabilidad
matriz_controlabilidad = [matriz_Bc matriz_Ac*matriz_Bc matriz_Ac^2*matriz_Bc matriz_Ac^3*matriz_Bc]; 
controlabilidad = rank(matriz_controlabilidad); 

num_variables_estado = size(estado_inicial);

if(num_variables_estado(1) == controlabilidad)
    disp('El sistema es controlable')
else
    disp('El sistema NO es controlable')
end

%%
% Etapa IV: Definición de la Ley de Control

%-------------------------------------------------------------------------%
%---------------------- Ley de Control con Integrador --------------------%
%-------------------------------------------------------------------------%
% Se establece una ley de control con integrador para anular el error en
% estado estable ante una referencia no nula.

%-------------------------------------------------------------------------%
%------------------------- Sistema Ampliado ------------------------------%
%-------------------------------------------------------------------------%
matriz_A_ampliada = [matriz_Ac  zeros(4,1) ; -matriz_Cc(1,:) 0];
matriz_B_ampliada = [matriz_Bc; 0];

%-------------------------------------------------------------------------%
%---------- Verificación de Controlabilidad del Sistema Ampliado ---------%
%-------------------------------------------------------------------------%
matriz_controlabilidad_ampliada = [matriz_B_ampliada matriz_A_ampliada*matriz_B_ampliada matriz_A_ampliada^2*matriz_B_ampliada matriz_A_ampliada^3*matriz_B_ampliada matriz_A_ampliada^4*matriz_B_ampliada]; 

num_variables_estado_ampliado = num_variables_estado(1)+1;

if(num_variables_estado_ampliado == rank(matriz_controlabilidad_ampliada))
    disp('El sistema ampliado es controlable')
else
    disp('El sistema ampliado NO es controlable')
end

%%
% Etapa V: Definición de Polos Deseados

%-------------------------------------------------------------------------%
%-------------------- Ubicación de Polos de Lazo Cerrado -----------------%
%-------------------------------------------------------------------------%

% Polos del sistema original
polos_original = [-15+15i, -15-15i, -0.5+0.5i, -0.5-0.5i];
% Polo del integrador
polo_integrador = -.000001;     

polos_deseados = [polos_original, polo_integrador];

%%
% Etapa VI: Implementación del Controlador

%-------------------------------------------------------------------------%
%------------------ Cálculo de la Matriz del Controlador -----------------%
%-------------------------------------------------------------------------%

% Transformación a Forma Controlable Canonica
valores_autos = eig(matriz_A_ampliada); 
coeficientes_autos = poly(valores_autos);

matriz_W = [coeficientes_autos(5) coeficientes_autos(4) coeficientes_autos(3) coeficientes_autos(2)  1;   
            coeficientes_autos(4) coeficientes_autos(3) coeficientes_autos(2)    1     0;   
            coeficientes_autos(3) coeficientes_autos(2)    1       0     0;
            coeficientes_autos(2)    1       0       0     0;
            1          0       0       0     0];

matriz_T = matriz_controlabilidad_ampliada * matriz_W;

matriz_A_controlable = inv(matriz_T) * matriz_A_ampliada * matriz_T;

coeficientes_polos_deseados = poly(polos_deseados);

matriz_K_fcc = fliplr([coeficientes_polos_deseados(2:end) - coeficientes_autos(2:end)]) * inv(matriz_T);

K_controlable_canonico = matriz_K_fcc(1:4);
K_integral_canonico = -matriz_K_fcc(end);

% Método de Ackerman
K_ackerman = acker(matriz_A_ampliada, matriz_B_ampliada, polos_deseados);
K_ack = K_ackerman(1:4);
K_integral_ack = -K_ackerman(end);

% LQR
pesos_q = [1 1000000 1 1 .1];     
matriz_Q = diag(pesos_q);
peso_R = 1000000;

K_lqr_total = lqr(matriz_A_ampliada, matriz_B_ampliada, matriz_Q, peso_R);
K_lqr = K_lqr_total(1:4);
K_integral_lqr = -K_lqr_total(5);

eig(matriz_A_ampliada - matriz_B_ampliada*K_lqr_total);

%%
% Etapa VII: Simulación

%-------------------------------------------------------------------------%
%---------------------------- Inicializaciones ---------------------------%
%-------------------------------------------------------------------------%

% Referencias
referencia = -100;

% Tiempos
tiempo_simulacion = 70;      % [s]
paso_integracion = 1e-4;    % [s]
num_pasos = tiempo_simulacion / paso_integracion;

% Vectores de Simulación
tiempo = 0:paso_integracion:(tiempo_simulacion-paso_integracion);
accion_control(1)  = 0;        
accion_control_efectiva(1) = 0;        

% Punto de Operación
punto_operacion =[0 0 0 referencia]';

% Integrador
integrador_error(1) = 0;

%%
%-------------------------------------------------------------------------%
%------------------------------- Simulación ------------------------------%
%-------------------------------------------------------------------------%

for i=1:num_pasos
   
   % Estados actuales
   estado_actual = [angulo_alpha(i); angulo_phi(i); vel_phi(i); altura_h(i)];
   
   % Salida actual
   salida_actual = matriz_Cc * estado_actual;
   
   % Error de Control
   error_control = referencia - salida_actual(1); 
   
   % Integral del Error de Control
   integrador_error(i+1) = integrador_error(i) + error_control * paso_integracion;
   
   % Acción de Control
   accion_control(i) = -K_lqr * estado_actual + K_integral_lqr * integrador_error(i+1);          
   
   % Sistema Lineal
   % Actualización de Variables
   derivada_estado = matriz_Ac * estado_actual + matriz_Bc * accion_control(i);      
   estado_actual = estado_actual + derivada_estado * paso_integracion;          
   
   if(i+1 <= num_pasos)
       angulo_alpha(i+1) = estado_actual(1);
       angulo_phi(i+1)   = estado_actual(2);
       vel_phi(i+1)      = estado_actual(3);
       altura_h(i+1)     = estado_actual(4);    
   end
  
end    

%%
%-------------------------------------------------------------------------%
%-------------------------------- Gráficas -------------------------------%
%-------------------------------------------------------------------------%

color_linea = 'b';

figure(1);

subplot(3,2,1); grid on; hold on;
plot(tiempo, angulo_alpha, color_linea);
title('Ángulo con la Horizontal, \alpha');
ylabel('\alpha [rad]');
xlabel('Tiempo [s]');

subplot(3,2,2); grid on; hold on;
plot(tiempo, angulo_phi, color_linea);
title('Ángulo de Cabeceo, \phi');
ylabel('\phi [rad]');
xlabel('Tiempo [s]');

subplot(3,2,3); grid on; hold on;
plot(tiempo, vel_phi, color_linea);
