clear 
clc
load('ecg118.mat');

A=r(1:6000,:);
x=A(:,2);

[F,C]=size(A);
y=zeros(F,1);

% Señal con filtro paso bajo
for i=13:F
   y(i)=2*y(i-1)-y(i-2)+(1/32)*[x(i)-2*x(i-6)+x(i-12)];
end

% Señal con filtro paso alto
for i=33:F
    y1(i)=x(i-1)-(y(i)/32)+y(i-16)-y(i-17)+(y(i-32)/32);
end

% Señal con derivador
for i=5:F
    y2(i)=(1/8)*[2*y1(i)+y1(i-1)-y1(i-3)-2*y1(i-4)];
end

% Señal con operación cuadratura
y3=y2.^2;

figure (1)
subplot(2,2,1)
plot(x)
title('Señal Original')
subplot(2,2,2)
plot(y)
title('Señal con filtro paso bajo')
subplot(2,2,3)
plot(y1)
title('Señal con filtro paso alto')
subplot(2,2,4)
plot(y2)
title('Señal con derivador')

% Señal con integrador de ventana móvil
N=30; %Valor adecuado para una frecuancia de muestreo de 200Hz
for i=31:length(y2)
    resto=0;
    for j=N:-1:0
       resto=resto+y3(i-j);
    end
    yI(i)=(1/N)*resto;
end

% Señal con umbral adaptativo
% picos SPK son picos de complejo QRS
% picos NPK son picos de ruido

[pks,locs]=findpeaks(yI);
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


for i=1:length(yI)
    if yI(i)<Umbral
        yf1(i)=0;
    else
        yf1(i)=yI(i);
    end
end

pasos=length(yf1)/300;
for i=1:pasos
    lims=300*i;
    limi=lims-299;
    C=yf1(limi:lims);
    c(i,2)=max(C);
    c(i,1)=find(yf1==c(i,2));
end

yf2=zeros(1,length(yf1));
for i=1:length(c)
    yf2(c(i,1))=1;
end

figure (2)
subplot(2,2,1)
plot(y3)
title('Señal con operación cuadratura')
subplot(2,2,2)
plot(yI)
title('Señal con integrador de ventana móvil')
subplot(2,2,3)
plot(yf1)
title('Señal filtro adaptativo')
subplot(2,2,4)
plot(yf2)
title('Muestras')

figure (3)
subplot(2,1,1)
plot(yI,'Color','b')
hold on
plot(yf2*2000,'Color','r');
title('Señal con integrador de ventana móvil y con complejos QRS identificados')
subplot(2,1,2)
plot(yf1,'Color','b')
for i=1:pasos
    viscircles([c(i,1) c(i,2)], 40,'Color','r');
end
title('Señal con filtro adaptativo y con complejos QRS identificados')



