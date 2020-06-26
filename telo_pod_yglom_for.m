% Движение под углом к горизонту с циклом for
clear all;

%начальная координата
x0=0; 
y0=0;

%начальная скорость
alpha=pi/4;
v0=150; 
vx0=v0*cos(alpha);
vy0=v0*sin(alpha);

k=30; %коэфф сопр
massa=30;
g=9.801;

%равномерная сетка по времени
tau=0.1; 
T=12;
M=T/tau;

%начальные условия
U(1,1)=x0;
U(1,2)=y0;
U(1,3)=vx0;
U(1,4)=vy0;

method = 1; %ERK_ 1,2,3,4

if method == 1
    for m=2:M
        w(m-1,1)=U(m-1,3);
        w(m-1,2)=U(m-1,4);
        w(m-1,3)=-k/massa*U(m-1,3);
        w(m-1,4)=-g - k/massa*U(m-1,4);
        for i=1:4
            U(m,i)=U(m-1,i)+tau*w(m-1,i);
        end
    end
end

if method==2 %ERK2
    a2=2/3;
    b1=1/4;
    b2=3/4;
    c1=0;
    c2=0;
    for m=2:M
        w1(m-1,1)=U(m-1,3);
        w1(m-1,2)=U(m-1,4);
        w1(m-1,3)=-k/massa*U(m-1,3);
        w1(m-1,4)=-g - k/massa*U(m-1,4);
        
        Uw1=U+tau*a2*w1;
        
        w2(m-1,1)=Uw1(m-1,3);
        w2(m-1,2)=Uw1(m-1,4);
        w2(m-1,3)=-k/massa*Uw1(m-1,3);
        w2(m-1,4)=-g - k/massa*Uw1(m-1,4);
        for i=1:4
            U(m,i)=U(m-1,i)+tau*(w1(m-1,i)*b1+b2*w2(m-1,i));
        end
    end
end

if method==3 %ERK3
    a2=1/2;
    a3=3/4; 
    b1=2/9;
    b2=3/9;
    b3=4/9;  
    for m=2:M
        w1(m-1,1)=U(m-1,3);
        w1(m-1,2)=U(m-1,4);
        w1(m-1,3)=-k/massa*U(m-1,3);
        w1(m-1,4)=-g - k/massa*U(m-1,4);
        
        Uw2=U+tau*a2*w1;
        
        w2(m-1,1)=Uw2(m-1,3);
        w2(m-1,2)=Uw2(m-1,4);
        w2(m-1,3)=-k/massa*Uw2(m-1,3);
        w2(m-1,4)=-g - k/massa*Uw2(m-1,4);
        
        Uw3=U+tau*a3*w2;
        
        w3(m-1,1)=Uw3(m-1,3);
        w3(m-1,2)=Uw3(m-1,4);
        w3(m-1,3)=-k/massa*Uw3(m-1,3);
        w3(m-1,4)=-g - k/massa*Uw3(m-1,4);
        for i=1:4
            U(m,i)=U(m-1,i)+tau*(w1(m-1,i)*b1+b2*w2(m-1,i)+b3*w3(m-1,i));
        end
    end
end

if method==4 %ERK4
    a2=1/2;
    a3=1/2; 
    a4=1;
    b1=1/6;
    b2=1/3;
    b3=1/3; 
    b4=1/6;
    for m=2:M
        w1(m-1,1)=U(m-1,3);
        w1(m-1,2)=U(m-1,4);
        w1(m-1,3)=-k/massa*U(m-1,3);
        w1(m-1,4)=-g - k/massa*U(m-1,4);
        
        Uw2=U+tau*a2*w1;
        
        w2(m-1,1)=Uw2(m-1,3);
        w2(m-1,2)=Uw2(m-1,4);
        w2(m-1,3)=-k/massa*Uw2(m-1,3);
        w2(m-1,4)=-g - k/massa*Uw2(m-1,4);
        
        Uw3=U+tau*a3*w2;
        
        w3(m-1,1)=Uw3(m-1,3);
        w3(m-1,2)=Uw3(m-1,4);
        w3(m-1,3)=-k/massa*Uw3(m-1,3);
        w3(m-1,4)=-g - k/massa*Uw3(m-1,4);
        
        Uw4=U+tau*a4*w3;
        
        w4(m-1,1)=Uw4(m-1,3);
        w4(m-1,2)=Uw4(m-1,4);
        w4(m-1,3)=-k/massa*Uw4(m-1,3);
        w4(m-1,4)=-g - k/massa*Uw4(m-1,4);
        for i=1:4
            U(m,i)=U(m-1,i)+tau*(w1(m-1,i)*b1+b2*w2(m-1,i)+b3*w3(m-1,i)+b4*w4(m-1,i));
        end
    end
end

figure(1);
plot(U(:,1),U(:,2),'.');
hold on
plot(U(:,1),0,'k-');
xlabel('X');
ylabel('Y');
title('Движение тела под углом к горизонту');
legend('Тело','Земля');