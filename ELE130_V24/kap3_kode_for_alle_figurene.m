%% ----------------------------------------------------------------
% figur 3.1
%----------------------------------------------------------------

clear; close all; clc

x = 0:pi/6:2*pi;
f = sin(x);

plot(x,f,'k*-')
xticks(x)
xticklabels({'0','\pi/6','2\pi/6','3\pi/6','4\pi/6','5\pi/6',...
    '\pi','7\pi/6','8\pi/6','9\pi/6','10\pi/6','11\pi/6','2\pi'})
grid
axis([0 2*pi -1.1 1.1])
xlabel('$x$ [rad]','FontSize',16)
ylabel('$f(x)$','FontSize',16)


%% ----------------------------------------------------------------
% figur 3.2
%----------------------------------------------------------------
clear; close all; clc

x = 0:0.01:2*pi;
a = 3;
c = pi/4;
d = 4;
k = 2;
f = a*sin(k*(x-c))+d;

figure
set(gcf,'Position',[1000 600 560 350])
plot(x,f,'b-')
xticks([0 pi 2*pi])
xticklabels({'0','\pi','2\pi'})

grid
axis([0 2*pi 0 8])
xlabel('$x$ [rad]','FontSize',16)
ylabel('$f(x)$','FontSize',16)
title('$f(x)=3\sin(2(x-\frac{\pi}{4}))+4$','FontSize',14)


%% ----------------------------------------------------------------
% figur 3.3
%----------------------------------------------------------------
clear; close all; clc

t = 0:pi/6:2*pi;
f = sin(t);

figure
set(gcf,'Position',[1000 600 560 290])
plot(t,f,'k*-')
xticks(round(t,2))
grid
axis([0 2*pi -1.1 1.1])
xlabel('tiden $t$ [s]','FontSize',16)
ylabel('$f(t)$','FontSize',16)


%% ----------------------------------------------------------------
% figur 3.4
%----------------------------------------------------------------
clear; close all; clc

t = 0:0.1:7;
f = sin(t);

figure
set(gcf,'Position',[1000 600 560 290])
plot(t,f,'b-')
hold on

t1 = 0:1:7;
f1 = sin(t1);
plot(t1,f1,'b*')

grid
axis([0 7 -1.1 1.1])
xlabel('tiden $t$ [s]','FontSize',16)
ylabel('$f(t)$','FontSize',16)


%% ----------------------------------------------------------------
% figur 3.5, versjon 1
%----------------------------------------------------------------
clear; close all; clc

t = 0:0.01:2*pi;
a = 3;
c = pi/4;
d = 4;
k = 2;
u = a*sin(k*t)+d;
y = a*sin(k*t- k*c)+d;

figure
set(gcf,'Position',[1000 600 560 350])
plot(t,u,'r-')
hold on
plot(t,y,'b-')

grid
axis([0 2*pi 0 8])
xlabel('tiden $t$ [s]','FontSize',16)
legend('$u(t)=3\sin(2t)+4$',...
    '$y(t)=3\sin(2t-\frac{\pi}{2})+4$','interpreter','latex')


%% ----------------------------------------------------------------
% figur 3.5, versjon 2
%----------------------------------------------------------------
clear; close all; clc

t = 0:0.01:2*pi;
u = 3*sin(2*t)+4;
y = 3*sin(2*t-pi/2)+4;

figure
set(gcf,'Position',[1000 600 560 350])
plot(t,y,'b-',t,u,'r')
grid
axis([0 2*pi 0 8])
plot(t,u,'r-')
hold on
plot(t,y,'b-')
grid
axis([0 2*pi 0 8])
xlabel('tiden $t$ [s]','FontSize',16)
legend('$u(t)=3\sin(2t)+4$',...
    '$y(t)=3\sin(2t-\frac{\pi}{2})+4$', ...
    'Location','best')


%% ----------------------------------------------------------------
% figur 3.6
%----------------------------------------------------------------
close all; clear; clc

delta_t = 0.01;
t_slutt = 1;
t=0:delta_t:t_slutt;
u=2*sin(5*t);

% utdrag for fargelegging
t2=0.4:delta_t:0.6;
u2=2*sin(5*t2);

figure
set(gcf,'Position',[302 525 1078 349])

subplot(1,2,1)
plot(t,u,'k')
hold on
area(t,u,'FaceColor','[0.7 0.7 0.7]')
area(t2,u2,'FaceColor','[0.5 0.5 0.5]')
% en svart strek i 0
plot([0 t(end)],[0 0],'k','LineWidth',2)
grid
title('$u(t)=2{\cdot}\sin(5t)$')
xlabel('tid [s]')
axis([0 1 -2.1 2.1])

subplot(1,2,2)
y=-0.4*cos(5*t)+0.4;
xlabel('tid [s]')
plot(t,y,'k-')
grid
title('$y(t)=-0.4{\cdot}\cos(5t)+0.4$')
xlabel('tid [s]')
axis([0 1 -1 1])


