%Algoritmo Pan-Tompkins
clear
close all
clc
% ecg100 2 muestras de 650000(canales) %plot(A(:,2)) %%Todas las filas de la columna 2
%Reconcortemos la señal, desde la 1 hasta la 6000 de la columna 2
load('ecg100.mat') 
x=A(1:6000,2) ;
subplot(2,3,1)
plot(x)
title('Señal Original')

%Filtro Pasabanda 
   s=6000;
 %Filtro Pasabajo
   x1=zeros(s,1);
   for n= 13:s %pues no puede haber un valor en cero
       x1(n)=2*x1(n-1)-x1(n-2)+(1/32)*[x(n)-2*x(n-6)+x(n-12)];
   end
   subplot(2,3,2)
   plot(x1)
   title('Señal Filtro Pasabajo')

    
   %Filtro PasaAto
    x2=zeros(s,1);
    for n= 33:s %pues no puede haber un valor en cero
         x2(n)=x2(n-1)-(x1(n)/32)+x1(n-16)-x1(n-17)+(x1(n-32)/32);
    end
   subplot(2,3,3)
   plot(x2)
   title('Señal Filtro Pasa Alto') 
   
%Derivador
x3=zeros(s,1);
for n=5:s
    x3(n)=(1/8)*[2*x2(n)+x2(n-1)-x2(n-3)-2*x2(n-4)];
end
subplot(2,3,4)
plot(x3)
title('Señal Derivada') 


%Operacion Cuadratica
x4=zeros(s,1);
for n=1:s
    x4(n)=x3(n).^2;
end
subplot(2,3,5)
plot(x4)
title('Señal Operacion Cuadratura')

%Para una frecuencia de muestreo de 200Hz, vamos a sumar (integrar) de a 30 muestras
%Integrador de ventana Movil
x5=zeros(s,1);
N=30;
for i=31:length(x3)
    k=0;
    for j=N:-1:0
        k=k+(x4(i-j));
    end
    x5(i)=(1/N)*k;
end
subplot(2,3,6)
plot(x5)
title('Señal Integrador de ventana Movil')


% Señal con umbral adaptativo
% picos SPK son picos de complejo QRS
% picos NPK son picos de ruido

[pks,locs]=findpeaks(x5);
B=sort(pks,'descend');
limite=B(1)/5;

for i=1:length(pks)
    if limite>B(i)
    pos=i;
    break
    end
end

SPK=mean(B(1:pos-1));
NPK=mean(B(pos:length(B)));
Umbral=(SPK+NPK)/2;

for i=1:length(pks)
    if pks(i)>Umbral
%         este pico es una señal de QRS
        SPK=0.125*pks(i)+0.875*SPK;
    else 
%         este pico es una señal de ruido
        NPK=0.125*pks(i)+0.875*NPK;
    end
Umbral=NPK+0.25*(SPK-NPK);
end


for i=1:length(x5)
    if x5(i)<Umbral
        x6(i)=0;
    else
        x6(i)=x5(i);
    end
end

pasos=length(x6)/300;
for i=1:pasos
    lims=300*i;
    limi=lims-299;
    C=x6(limi:lims);
    c(i,2)=max(C);
    c(i,1)=find(x6==c(i,2));
end

yf2=zeros(1,length(x6));
for i=1:length(c)
    yf2(c(i,1))=1;
end

figure
plot(x6,'Color','b')
for i=1:pasos
    viscircles([c(i,1) c(i,2)], 40,'Color','r');
end
title('Señal con Umbral Adaptativo')


figure
imshow(imcomplement(x6))

