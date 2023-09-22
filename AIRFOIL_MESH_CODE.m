%LIMITES
a=0.0; % Inicial %(En a se va a mallar el perfil aerodinámico)
b=0.3+a ; % FInal (La capa exterior de la capa)
v0=0 ; %Punto inicial para que de la vuelta
v1=2*pi ; %Punto Final para que de la vuelta
 
% NUMERO DE NODOS
M= 30 ; %Numero de curvas que son paralelos al perfil aerodinámico
N = 1000 ; %Número de curvas que son ortogonales alperfil aerodinámico

%PARAMETROS
e0=0.75;
 %e0 := (Distancia del centro del círculo al origen)/(Distancia del punto donde el mapeo de Joukowski no es regular
 %SUGERENCIA: no hacer e0 muy grande
phi=368*pi/1000 ; %phi := ángulo con el que se aleja el círculo al origen de coordenadas
% IMPORTANTE pi/2 <= phi <= pi


%COORDENADAS Mapeamos (x,y) ---> (d, v).
%d:= la coordenada que se aleja del perfil aerodinámico
%v:= la coordenada que recorre el perfil aerodinámico

%Partición para los puntos a lo largo de d (Exponenicial)
%Modular la partición exponencial c
for i=1 : M
    d(i) = a + (b-a)*part((i-1)/(M-1)) ;
end

%Partición para los puntos a lo largo de v
for j=1 : N
    v(j) = v0 + (v1-v0)*cospart((j-1)/(N-1)) ;    %cospart((j-1)/(N-1)) ;
end


%MAPEO: Se definen x,y en términos de d,v
for j=1 : N
for i=1 : M
x(i,j) = Af(d(i), v(j),e0,phi) ;
y(i,j) = Ort(d(i), v(j),e0,phi) ;
end
end



%CUERDA (HIPOTESIS)
N2 = 100 ;
for i=1 : N2
   l(i) = -2 + 2*part((i-1)/(N2-1)) ;
xl(i) = Af(l(i), 0,e0,phi) ;
yl(i) = Ort(l(i), 0,e0,phi) ;
end

%PARA TORCER LA GEOMETRÍA UN ÁNGULO THETA (EN RADIANES)
%theta = pi/3

%for j=1 : N
%for i=1 : M
%x(i,j) = Af(d(i), v(j),e0,phi)*cos(theta) -  ;
%y(i,j) = Ort(d(i), v(j),e0,phi) ;
%end
%end



%GRAFICA (Falta darle formato)
hold on
%plot (x(1,:), y(1,:), 'LineWidth',2,'Color',[0,0,0]);
%for i=2 :M
%plot (x(i,:), y(i,:), 'LineWidth',1.5,'Color',[0,0,1]);
%end
%for j=1 :N
%plot (x(:,j), y(:,j), 'LineWidth',1.5,'Color',[1,0,0]) ;
%end

for k=1 : N2
plot (xl(k), yl(k), 'LineWidth',1.0,'Color',[0,0,0], 'LineStyle','-') ;
end

title('Curva con \epsilon = 0.25 y \phi_{0} = -\pi/2. ' ) 
ax = gca;
ax.TitleFontSizeMultiplier = 5;
xlim([-4,4])
ylim([-4,4])
hold off


G=[x(1,:); y(1,:)];
G=transpose(G) ;

%Par de funciones que generan las coordenadas:
% x = x(d,v;e0;phi)
% y = y(d,ve0;phi)
% e0 y phi son parametros. 
% e0 define el grosor del perfil 
% phi define la curvatura del perfil
function x=Af(d,v,e0,phi) ;
x=1.2*( (1+d)*(cos(v) +e0*cos(v - phi)) - e0*cos(phi) )*(1+ 1/( (1+d)^2*(1+e0^2+2*e0*cos(phi)) -2*(1+d)*e0*(cos(v+phi)+e0*cos(v) ) + e0^2 )) ;
% (1+ 1/( (1+d)^2 +2*(1+d)*e0*(cos(v-phi)-cos(phi) ) + 2*e0^2*(1-cos(v)) ))
%( (1+d)*cos(v) + e0*(cos(phi) - cos(v + phi)) )
end

function y=Ort(d,v,e0,phi) ;
y=1.2*( (1+d)*(sin(v) +e0*sin(v - phi)) + e0*sin(phi) )*(1 - 1/( (1+d)^2*(1+e0^2+2*e0*cos(phi)) -2*(1+d)*e0*(cos(v+phi)+e0*cos(v) ) + e0^2 )) ;
end


function w=part(i)
w =(1.2^i-1)/(1.2-1);
end

function z=cospart(i)
z = (1-cos(pi*i))/2 ;
end