%-----------------------------------------------------------------------
%% a) foregår i Simulink

%%

%-----------------------------------------------------------------------
%% b) Numerisk integrasjon og derivasjon i Matlab
% Sammenligning mot analytiske uttrykk
clear; close all; clc


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Starter først med "fasiten", hvor vi definerer kontinuerlig "tid"
% og analytiske uttrykk for u(t), y(t) og v(t) tilsvarende
% ligningene (3), (4) og (5) i oppgaven. Det vil si:
%
% u(t) = 2*t^2
% y(t) = 2/3*t^3, integrert av u(t)
% v(t) = 4*t, derivert av u(t)
%
% Liten verdi på "delta_t" her gir tilnærmet kontinuerlig tid.
% Benytter syntaks "tid" for å skille det fra diskret vektor "t"
% som brukes for t(k).
delta_t = 0.0001;
t_slutt = 3;
tid = 0:delta_t:t_slutt;

% Benytter syntaks "u_t", "y_t" og "v_t" for å fremheve at de
% analytiske uttrykkene er tidskontinuerlige.
% Husk elementvise operasjoner!!!!!!!!

u_t = 2*tid.^2;         % analytisk uttrykk for u(t)
y_t = 2/3*tid.^3;         % analytisk uttrykk for y(t)
v_t = 4*tid;         % analytisk uttrykk for v(t)


% Plotter for å sjekke at det ser riktig ut
figure
subplot(3,1,1)
plot(tid, y_t, 'b-')    % "fasit" analytisk uttrykk
grid on
title('Signal $y(t)$ og numerisk integrert $y_{k}$')

subplot(3,1,2)
plot(tid, u_t, 'k-')
grid on
title('Signal $u(t)$ og målinger $u_{k}$')

subplot(3,1,3)
plot(tid, v_t, 'r-')      % "fasit" analytisk uttrykk
grid on
xlabel('tid [s]')
title('Signal $v(t)$ og numerisk derivert $v_{k}$')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% Diskret tids- og signalvektor. Plukker ut målinger fra signalet u(t)
% ved å definere en tidsvektor t = [0, Ts, 2*Ts, ..., t_slutt]
% Benytter syntaks "t" og "u" for de diskret vektorene $t_k$ og $u_k$.
% Indeksering av t og u blir da som t(k) og u(k)
T_s = 0.4;              % samme Ts som i oppgave a) i Simulink
t_slutt = 3;
t = 0:T_s:t_slutt;      % diskret tidsvektor $t_k$
u = 2*t.^2;             % diskret målinger $u_k$


% Initialverdier integrasjon
y_EulerF(1) = 0;
y_EulerB(1) = 0;
y_Trapes(1) = 0;

% Initialverdier derivasjon
v_bakover(1) = 0;
v_senter(1) = 0;

for k = 2:length(u)
    % - - - - Integrasjonsrutinene - - - -
    y_EulerF(k) = y_EulerF(k-1) + T_s * u(k-1);
    y_EulerB(k) = y_EulerB(k-1) + T_s * u(k);
    y_Trapes(k) = y_Trapes(k-1) + T_s * 0.5 *(u(k-1)+u(k));

    % - - - - Derivasjonsrutinene - - - -
    v_bakover(k) = (u(k)-u(k-1)) / T_s;
    v_forover(k-1) = (u(k)-u(k-1)) / T_s;
    if k > 2
        v_senter(k-1) = (u(k)-u(k-2)) / (2*T_s);
    end
end

% Går tilbake til figuren, starter med hold on
subplot(3,1,1)
hold on
plot(t(1:k), y_EulerF, 'b:o')
plot(t(1:k), y_EulerB, 'b:v')
plot(t(1:k), y_Trapes, 'b:.','MarkerSize',15)
legend('Analytisk','Eulers forovermetode', ...
    'Eulers bakovermetode','Trapesmetoden',...
    'Location','northwest')

subplot(3,1,2)
hold on
plot(t, u, 'ko')

subplot(3,1,3)
hold on
plot(t(1:k), v_bakover, 'r:o')
plot(t(1:k-1), v_forover, 'r:v')
plot(t(1:k-1), v_senter, 'r:.','MarkerSize',15)

% Calculate the index k for t(k)=2
k = round(2 / T_s) + 1;

% Plot the value of v_bakover(k) at t(k)=2
subplot(3,1,3)
plot(t(k), v_bakover(k), 'go', 'MarkerSize', 10)

% Annotate the value of v_bakover at t(k)=2
text(t(k), v_bakover(k), [' (' num2str(t(k)) ', ' num2str(v_bakover(k)) ')'])

axis('tight')
legend('Analytisk', 'Bakoverderivasjon',...
    'Foroverderivasjon',...
    'Senterderivasjon',...
    'v_bakover(2)',...
    'Location','northwest')

%{
%-----------------------------------------------------------------------
%% c) foregår for hånd





%-----------------------------------------------------------------------
%% d) Numerisk integrasjon som Matlab-funksjon
clear;close all;clc

% Diskret tids- og signalvektor t(k) og u(k).
% Kopier tilsvarende kode fra deloppgave b).
T_s = ..;
t_slutt = 3;
t = ..;
u = ..;

% initialverdi
y(1) = ..;
for ..
    y(k) = EulersForover(..);
end

figure
subplot(2,1,1)
plot(t, u, 'k:o')
grid
title('Signal $u_{k}$')

subplot(2,1,2)
plot(t, y, 'b:o')
grid
title('$y_{k}$, numerisk integrert av $u_{k}$')
legend('Eulers forovermetode')


%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% d), generell funksjon for numerisk integrasjon (frivillig)
clear;close all;clc

% Diskret tids- og signalvektor t(k) og u(k).
% Kopier tilsvarende kode fra deloppgave b).
T_s = ..;
t_slutt = 3;
t = ..;
u = ..;

% initialverdi
y1(1) = ..;
y2(1) = ..;
y3(1) = ..;

for ..
    y1(k) = Integrasjon(y1(k-1), T_s, .. , metode='EulersForover');
    y2(k) = Integrasjon(y2(k-1), T_s, .. , metode='EulersBakover');
    y3(k) = Integrasjon(y3(k-1), T_s, .. , metode='Trapes');
end

figure
set(gcf,"Position",[180 130 600 800])

subplot(4,1,1)
% siste element skal ikke være med
bar(t(1:end-1)+T_s/2,u(1:end-1),1,'FaceColor',[0.7 0.7 0.7])
hold on
plot(t,u,'k:o')
grid
title('Estimert arealet under $u_{k}$ ved Eulers forovermetode')
axis([0 3 -5 22])

subplot(4,1,2)
bar(t-T_s/2,u,1,'FaceColor',[0.7 0.7 0.7])
hold on
plot(t,u,'k:o')
grid
title('Estimert arealet under $u_{k}$ ved Euler bakovermetode')
axis([0 3 -5 22])

subplot(4,1,3)
area(t,u,'FaceColor',[0.7 0.7 0.7])
hold on
stem(t,u,'k')
grid
title('Estimert arealet under $u_{k}$ ved trapesmetoden')
axis([0 3 -5 22])

subplot(4,1,4)
plot(t, y1, 'b:o')
hold on
plot(t, y2, 'b:v')
plot(t, y3, 'b:.','markersize',15)
grid
title('$y_{k}$, numerisk integrert av $u_{k}$')
legend('Eulers forovermetode', 'Eulers bakovermetode',...
    'Trapesmetoden','location','northwest')





%-----------------------------------------------------------------------
%% e) Numerisk derivasjon som Matlab-funksjon
clear;close all;clc

% Diskret tids- og signalvektor t(k) og u(k).
% Kopier tilsvarende kode fra deloppgave b).
T_s = ..;
t_slutt = 3;
t = ..;
u = ..;

%initialverdi
v(1) = ..;
for ..
    v(k) = BakoverDerivasjon(.. , ..);
end

figure
subplot(2,1,1)
plot(t, u,'k:o')
grid
title('Signal $u_{k}$')

subplot(2,1,2)
plot(t, v, 'r:o')
grid
title('$v_{k}$, numerisk derivert av $u_{k}$')
legend('Bakoverderivasjon')


%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% e) generell funksjon for numerisk derivasjon (frivillig)
clear; close all;clc

% Diskret tids- og signalvektor t(k) og u(k).
% Kopier tilsvarende kode fra deloppgave b).
T_s = ..;
t_slutt = 3;
t = ..;
u = ..;

% initialverdi
v1(..) = ..;   % Bakoverderivasjon
v2(..) = ..;     % Foroverderivasjon
v3(..) = ..;     % Senterderivasjon

for k = ..:length(t)
    v1(k)   = Derivasjon(u(k-2:k), T_s, metode='Bakover');
    v2(k-1) = Derivasjon(u(k-2:k), T_s, metode='Forover');
    v3(k-1) = Derivasjon(u(k-2:k), T_s, metode='Senter');
end

figure
subplot(2,1,1)
plot(t, u,'k:o')
grid
title('Signal $u_{k}$')

subplot(2,1,2)
plot(t, v1, 'r:o')
hold on
plot(t(1:end-1), v2, 'r:v')
plot(t(1:end-1), v3, 'r:.','markersize',15)
grid
title('$v_{k}$, numerisk derivert av $u_{k}$')
legend('Bakoverderivasjon','Foroverderivasjon',...
    'Senterderivasjon','location','southeast')





%-----------------------------------------------------------------------
%% f) Integrasjon av et forsterket signal
clear;close all;clc

% Bestemmer t(k) og u(k)
T_s = 0.4;
t_slutt = 3;
t = 0:T_s:t_slutt;
u = ones(1,length(t));
%u(5:end)=0;

% initialverdi
y1(1) = 0;
y2(1) = 0;

Ki = 0.4;
for k = ..
    y1(k) =    EulersForover(y1(k-1), T_s, Ki*u(k-1));
    y2(k) = Ki*EulersForover(y2(k-1), T_s,    u(k-1));
end

figure
set(gcf,'position',[200 70 600 700])
subplot(3,1,1)
bar(t+T_s/2,u,1,'FaceColor',[0.7 0.7 0.7])
hold on
plot(t, u,'ko')
grid
title('Signal $u_{k}$ og arealet under kurven ved Eulers forovermetode' )
yl = ylim;
ylim(yl*1.3)
xlim([0 3])

subplot(3,1,2)
bar(t+T_s/2,Ki*u,1,'FaceColor',[0.7 0.7 0.7])
hold on
plot(t, Ki*u,'ko')
grid
title(['Signal $K_i{\cdot}u_{k}$ og arealet under',...
    ' kurven ved Eulers forovermetode, $K_i$=',num2str(Ki)])
yl = ylim;
ylim(yl*1.3)
xlim([0 3])

subplot(3,1,3)
plot(t, y1, 'b:o')
hold on
plot(t, y2, 'r:o')
grid
title('Signalene $y_{1,k}$ og $y_{2,k}$')
yl = ylim;
ylim(yl*1.3)
legend('{\tt y1(k) = EulersForover(y1(k-1), T\_s, Ki*u(k-1))}',...
    '{\tt y2(k) = Ki*EulersForover(y2(k-1), T\_s, u(k-1))}',...
    'location','northwest')






%-----------------------------------------------------------------------
%% g) Effekten av støy ved numerisk derivasjon
clear;close all;clc

% Bestemmer t(k) og u(k)
T_s = 0.4;
t_slutt = 3;
t = 0:T_s:t_slutt;
u1 = ones(1,length(t));
u2 = u1 + 0.02*randn(1,length(t));

% initialverdi
y1(1) = 0;     % Eulers forovermetode på signal uten støy
y2(1) = 0;     % Eulers forovermetode på signal med støy

v1(..) = 0;   % bakoverderivasjon av signal uten støy
v2(..) = 0;   % bakoverderivasjon av signal med støy
v3(..) = 0;     % senterderivasjon av signal med støy
for k = 2:length(t)
    y1(k) = Integrasjon(y1(k-1), T_s, u1(k-1:k), metode='EulersForover');
    y2(k) = Integrasjon(y2(k-1), T_s, u2(k-1:k), metode='EulersForover');

    if k ..
        v1(k)   = Derivasjon(u1(k-2:k), T_s, metode='Bakover');
        v2(k)   = Derivasjon(u2(k-2:k), T_s, metode='Bakover');
        v3(k-1) = Derivasjon(u2(k-2:k), T_s, metode='Senter');
    end
end

figure
set(gcf,'position',[100 200 800 600])
subplot(3,2,1)
plot(t, u1, 'k:o')
grid
title('Signal {\tt u1}')
legend('Signal uten st{\o}y',...
    'location','northeast')
ylim([0.8 1.2])

subplot(3,2,2)
plot(t, u2,'k:o')
grid
title('Signal {\tt u2}')
legend('Signal med st{\o}y',...
    'location','northeast')
ylim([0.8 1.2])

subplot(3,2,3)
plot(t, v1, 'r:o')
grid
title('Signal {\tt v1}')
legend('Bakoverderivasjon av signal uten st{\o}y',...
    'location','northeast')
ylim([-0.5 0.5])

subplot(3,2,4)
plot(t, v2, 'r:o')
grid
title('Signal {\tt v2}')
legend('Bakoverderivasjon av signal med st{\o}y',...
    'location','northeast')
ylim([-0.1/T_s 0.1/T_s])

subplot(3,2,5)
plot(t, y1, 'b:o')
hold on
plot(t, y2, 'b:.')
xlabel('tiden $t$ [s]')
grid
title('Signal {\tt y1} og {\tt y2}')
legend('Eulers forovermetode av signal uten st{\o}y',...
    'Eulers forovermetode av signal med st{\o}y',...
    'location','best')

subplot(3,2,6)
plot(t(1:end-1), v3, 'r:.','markersize',15)
title('Signal {\tt v3}')
xlabel('tiden $t$ [s]')
grid
legend('Senterderivasjon av signal med st{\o}y',...
    'location','northeast')
ylim([-0.1/T_s 0.1/T_s])

sgtitle(['Tidsskritt $T_s$=',num2str(T_s),' sekund'])






%-----------------------------------------------------------------------
%% h) Fenomenet nedfolding
clear; close all; clc

f_1 = 0.1;  % Hz
f_2 = 0.3;  % Hz
f_3 = 0.5;  % Hz
f_4 = 0.8;  % Hz
f_5 = 1.3;  % Hz

% "kontinuerlig" tid for analoge signal
t_slutt = 15; % sekund
delta_t = 0.001;
tid = 0:delta_t:t_slutt;  % liten delta_t gir "kontinuerlig" tid

% Bestemmer bestemmer t(k) basert på samplingsfrekvens
f_s = 0.6;     % [Hz], samplingsfrekvens, endre på denne
T_s = ..;   % [sek], sampletid, beregnet fra f_s
t = 0:T_s:t_slutt;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Lager et "analogt" og et diskret signal av de fem signalene
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x1_a = sin(2*pi*f_1*tid);   % "analogt" signal
x1_d = sin(2*pi*f_1*t);     % diskret signal

x2_a = sin(2*pi*f_2*tid);
x2_d = sin(2*pi*f_2*t);

x3_a = sin(2*pi*f_3*tid);
x3_d = sin(2*pi*f_3*t);

x4_a = sin(2*pi*f_4*tid);
x4_d = sin(2*pi*f_4*t);

x5_a = sin(2*pi*f_5*tid);
x5_d = sin(2*pi*f_5*t);

% Lager et sammensatt signal av alle
% de 5 signalene
x_a = x1_a + x2_a + x3_a + x4_a + x5_a;
x_d = x1_d + x2_d + x3_d + x4_d + x5_d;

figure
set(gcf,'position',[550 100 550 850])
subplot(7,1,1)
plot(tid, x1_a,'r-', t, x1_d,'b:o')
title(['$x_1(t)=\sin(2\pi',num2str(f_1),'t)$'])

subplot(7,1,2)
plot(tid, x2_a,'r-', t, x2_d,'b:o')
title(['$x_2(t)=\sin(2\pi',num2str(f_2),'t)$'])

subplot(7,1,3)
plot(tid, x3_a,'r-', t, x3_d,'b:o')
title(['$x_3(t)=\sin(2\pi',num2str(f_3),'t)$'])

subplot(7,1,4)
plot(tid, x4_a,'r-', t, x4_d,'b:o')
title(['$x_4(t)=\sin(2\pi',num2str(f_4),'t)$'])

subplot(7,1,5)
plot(tid, x5_a,'r-', t, x5_d,'b:o')
title(['$x_5(t)=\sin(2\pi',num2str(f_5),'t)$'])

sgtitle(['Sampling av analoge signal. ',...
    'Samplingsfrekvens $f_s=$ ',num2str(f_s),' Hz'])

subplot(7,1,6)
plot(tid, x_a,'r-', t, x_d,'bo')
ax = axis;
title(['Analogt signal $x(t)=\sum_{i=1}^5 x_i(t)$ ',...
    'inneholder frekvensomr{\aa}det ',...
    num2str(f_1),'-',num2str(f_5),' Hz'])

subplot(7,1,7)
plot(t, x_d,'b-')
axis(ax)
title(['Med $f_s{=}$',num2str(f_s),' Hz er maksimal',...
    ' frekvens i $x_{k}$ lik $f_N{=}$',num2str(f_s/2),' Hz'])





%-----------------------------------------------------------------------
%% i) FIR-filter
clear; close all; clc

% Bestemmer t(k)
f_s = 10;
T_s = 1/f_s;
t_slutt = 10;
t = 0:T_s:t_slutt;

% Målesignalet x(k) er bare støy
rng(8)
x = rand(1,length(t));

% initialisering
y_FIR(1) = ..;

for k = ..
    M = 10;
    % juster på M her
    y_FIR(k) = ..
end

figure
plot(t, x, 'b','linewidth',1)
hold on
plot(t, y_FIR, 'r-','linewidth',1)
title(['FIR-filtrering, $T_s=$ ',num2str(T_s)])
legend('signal $x_k$',...
    ['$y_{k,FIR}$, $M=$ ',num2str(M)])


%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% i) FIR-filter, effekten av M
clear; close all; clc

% Bestemmer t(k)
f_s = 10;
T_s = 1/f_s;
t_slutt = 10;
t = 0:T_s:t_slutt;

% Målesignalet x(k) er bare støy
rng(2)
x = rand(1,length(t));

% initialisering
y_FIR1(1) = ..
y_FIR2(1) = ..
y_FIR3(1) = ..
M = [10, 30, 50];

for k = 2:length(t)
    y_FIR1(k) = FIR_filter(x(1:k),M(1));
    y_FIR2(k) = FIR_filter(x(1:k),M(2));
    y_FIR3(k) = FIR_filter(x(1:k),M(3));
end

figure
plot(t, x, 'b','linewidth',1)
hold on
plot(t, y_FIR1, 'r-','linewidth',1)
plot(t, y_FIR2, 'g-','linewidth',1)
plot(t, y_FIR3, 'k-','linewidth',1)
title(['FIR-filtrering, $T_s=$ ',num2str(T_s)])
legend('signal $x_k$',...
    ['$y_{k,FIR}$, $M=$ ',num2str(M(1))],...
    ['$y_{k,FIR}$, $M=$ ',num2str(M(2))],...
    ['$y_{k,FIR}$, $M=$ ',num2str(M(3))])




%-----------------------------------------------------------------------
%% j) IIR-filter
clear; close all; clc

% Bestemmer t(k)
f_s = 10;
T_s = 1/f_s;
t_slutt = 10;
t = 0:T_s:t_slutt;

% målesignalet x(k) er bare støy
rng(8)
x = rand(1,length(t));

% initialisering
y_IIR(1) = ..
alfa = 0.05;

for k = 2:length(t)
    y_IIR(k) = ..
end

figure
plot(t, x, 'b','linewidth',1)
hold on
plot(t, y_IIR, 'r-','linewidth',1)
title(['IIR-filtrering, $T_s=$ ',num2str(T_s)])
legend('signal $x_k$',...
    ['$y_{k,IIR}$, $\alpha=$ ',num2str(alfa)])


%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% j) IIR-filter, effekten av alfa
clear; close all; clc

% Bestemmer t(k)
f_s = 10;
T_s = 1/f_s;
t_slutt = 10;
t = 0:T_s:t_slutt;

% Målesignalet x(k) er bare støy
rng(8)
x = rand(1,length(t));

% initialisering
y_IIR1(1) = ..
y_IIR2(1) = ..
y_IIR3(1) = ..
alfa = [0.05, 0.1, 0.2];

for k = 2:length(t)
    y_IIR1(k) = IIR_filter(y_IIR1(k-1),x(k),alfa(1));
    y_IIR2(k) = IIR_filter(y_IIR2(k-1),x(k),alfa(2));
    y_IIR3(k) = IIR_filter(y_IIR3(k-1),x(k),alfa(3));
end


figure
plot(t, x, 'b','linewidth',1)
hold on
plot(t, y_IIR1, 'r-','linewidth',1)
plot(t, y_IIR2, 'g-','linewidth',1)
plot(t, y_IIR3, 'k-','linewidth',1)
title(['IIR-filtrering, $T_s=$ ',num2str(T_s)])
legend('signal $x_k$',...
    ['$y_{k,IIR}$, $\alpha=$ ',num2str(alfa(1))],...
    ['$y_{k,IIR}$, $\alpha=$ ',num2str(alfa(2))],...
    ['$y_{k,IIR}$, $\alpha=$ ',num2str(alfa(3))])



%-----------------------------------------------------------------------
%% k) Filtrering av ulike testsignal
clear; close all; clc

% Bestemmer t(k)
f_s = 2;      % samplingsfrekvens Hz
T_s = 1/f_s;  % samplingstid sekund
t_slutt = 30;
t = 0:T_s:t_slutt;

% - - - - - - - - - - - - - - - - - - - -
% x1, impuls
% - - - - - - - - - - - - - - - - - - - -
x1 = zeros(1,length(t));
t_impuls = 10;        % tidspunkt impuls
k = t_impuls/T_s+1;   % k*Ts=t_impuls, '+1' pga Matlabindeks
x1(k) = 1;            % en 1'er ved t=t_impuls sekund

% - - - - - - - - - - - - - - - - - - - -
% x2, firkantpuls (sprang opp og ned)
% - - - - - - - - - - - - - - - - - - - -
x2 = zeros(1,length(t));
t_sprang1 = 10;        % tidspunkt sprang start
t_sprang2 = 20;        % tidspunkt sprang slutt
k1 = t_sprang1/T_s+1;  % k1*Ts=t_sprang1
k2 = t_sprang2/T_s+1;  % k2*Ts=t_sprang2
x2(k1:k2) = 1;

% - - - - - - - - - - - - - - - - - - - -
% x3, sinus med frekvens f=0.1 Hz
% - - - - - - - - - - - - - - - - - - - -
f = 0.1;
x3 = 0.1*sin(2*pi*f*t) + 0.2;

% - - - - - - - - - - - - - - - - - - - -
% x4, sinus med støy
% - - - - - - - - - - - - - - - - - - - -
x4 = x3 + 0.03*randn(1,length(t));

% ------------------------------------------------
% Velg hvilket signal som skal filtreres.
x = x1;
% For x1 og x2 er det smart å bruke kurveform = ':x';
% For x3 og x4 er det smart å bruke kurveform = '-';
kurveform = ':x';
% ------------------------------------------------

% initialisering
y_FIR(1) = ..
y_IIR(1) = ..
alfa = 0.3;
M = 10;

for k = 2:length(t)
    y_FIR(k) = FIR_filter(x(1:k),M);
    y_IIR(k) = IIR_filter(y_IIR(k-1),x(k),alfa);
end

figure
subplot(2,2,1)
plot(t, x, ['b',kurveform],'linewidth',1)
grid
title(['signal $x_k$, $f_s=$ ',num2str(f_s),...
    ' Hz, $T_s=$ ',num2str(T_s),' s'])
yl = ylim;
yl = [yl(1)*0.9, yl(2)*1.05];
ylim(yl)

subplot(2,2,2)
plot(t, y_FIR, ['r',kurveform],'linewidth',1)
grid
title(['$y_{k,FIR}$, $M=$ ',num2str(M)])
ylim(yl)

subplot(2,2,3)
plot(t, y_IIR, ['g',kurveform],'linewidth',1)
grid
title(['$y_{k,IIR}$, $\alpha=$ ',num2str(alfa)])
ylim(yl)

subplot(2,2,4)
plot(t, x, ['b',kurveform],'linewidth',1)
hold on
plot(t, y_FIR, ['r',kurveform],'linewidth',1)
grid
title('Sammenligning')
ylim(yl)



%-----------------------------------------------------------------------
%% l) IIR-basert lavpass- og høypassfiltrering av Chirp-signal
clear; close all; clc

% Bestemmer t(k)
f_s = 1000;
T_s = 1/f_s;
t_slutt = 10;
t = 0:T_s:t_slutt;

f_initial = 0.1;
f_target = 10;
t_target = t_slutt;
f = ..
x = ..

y1_IIR(1) = x(1);
alfa = 0.01;

y2_IIR(1) = 0;
a1 = 0.95;
b0 = 0.7;
b1 = -0.7;

for k = 2:length(t)
    y1_IIR(k) = (1-alfa)*y1_IIR(k-1) + alfa*x(k);
    y2_IIR(k) = a1*y2_IIR(k-1) + b0*x(k) + b1*x(k-1);
end

figure
plot(t, x, 'b')
hold on
plot(t, y1_IIR, 'r')
plot(t, y2_IIR, 'g')
title(['IIR-filtrering av chirp-signal, $f_s=$ ',num2str(f_s),' Hz'])
legend('signal $x_k$',...
    ['Lavpass: $y_{k,IIR}$, $\alpha=$ ',num2str(alfa)],...
    ['H{\o}ypass: $y_{k,IIR}$, $a_1=$ ',num2str(a1),...
    ', $b_0=$ ',num2str(b0),', $b_1=$ ',num2str(b1)])
yl = ylim;
ylim(yl*1.1)


%-----------------------------------------------------------------------
%% m) foregår i Simulink




%-----------------------------------------------------------------------
%% n) Demonstrasjon av sampling og nedfolding
%
% Innhold:
%  - Visualisering av sampling i "sanntid"
%  - Velg frekvens f og samplingsfrekvens f_s fra meny
%
% ---------------------------------------------------------

clear; close all; clc

figure
set(gcf,'position',[400 400 650 700])
subplot(3,1,1)

% Høyeste frekvens [Hz] på signalet y(t) som skal samples.
% Juster gjerne på denne og studer effekten
% --------------------------------------------------------
f = 4.2;
% --------------------------------------------------------
T_p = 1/f;           % periode for signalet
delta_t = T_p/20;    % kontinuerlig tid består av 20 punkter i en periode
t_slutt = 6;
tid = 0:delta_t:t_slutt;  % kontinuerlig tid for "analogt" signal

%------------------------------------------------------------
% Velg hvilket signal du vil sample
% Du kan også lage helt egne signal om du vil.
%------------------------------------------------------------
y_t = sin(2*pi*f*tid);
%y_t = sin(2*pi*f*tid) + sin(2*pi*0.5*f*tid);
%y_t = sin(2*pi*f*tid) + sin(2*pi*0.5*f*tid) + sin(2*pi*0.1*f*tid);


% Plotter i sin helhet det analoge signalet yt) som skal samples
plot(tid,y_t,'r')
axis([0 t_slutt min(y_t)*1.1 max(y_t)*1.1])
title(['Analogt sinussignal $y(t)$ med høyeste frekvens $f=',...
    num2str(f),'$ Hz som skal samples.'])

% Initialverdi for samplet signal y(k)
y(1) = y_t(1);

% Indeks for sampling nummer k
k = 1;

% Åpner en dialogboks for input. Svaret fs er en cell
fs = inputdlg('Hvilken samplingsfrekvens vil du bruke? ');
Ts = 1/eval(fs{:});
% lager samplingstidspunkt t(k)
t_k = 0:Ts:t_slutt+Ts;

for i = 1:length(tid)

    % Viser det analoge signalet y(t) i "sanntid" frem til tid(i)
    subplot(3,1,2)
    plot(tid(1:i),y_t(1:i),'r')
    axis([0 t_slutt min(y_t)*1.1 max(y_t)*1.1])
    hold on
    grid on
    title(['Signal $y(t)$ og m{\aa}linger $y_k$,',...
        'samplingsfrekvens $f_s=',num2str(1/Ts), '$'])

    % Sjekker om vi skal sample
    if tid(i) >= t_k(k)
        % tar en måling y(k) av analogt signal y_t
        y(k) = y_t(i);
        k = k + 1;
    end
    plot(t_k(1:k-1),y(1:k-1),'bo')
    drawnow

    % benytter hold off for å unngå at figuren består av 100-vis av linjer
    hold off

    subplot(3,1,3)
    plot(t_k(1:k-1),y(1:k-1),'bo')
    axis([0 t_slutt min(y_t)*1.1 max(y_t)*1.1])
    grid on
    title('M{\aa}linger $y_k$')
    xlabel('tiden $t$ [s]')

end

hold on
plot(t_k(1:k-1),y(1:k-1),'b-')



% ---------------------------------------------------------
%% o)
%
% Innhold:
%  - Oppsamling ved bruk av interp-funksjonen
%  - beregning av sinus med sind-funksjonen
%
% ---------------------------------------------------------


clear; close all; clc

% Fakta om vinkelmålingen
delta_phi = 220;      % antall grader mellom hver måling
T_s = 0.2;            % tidsskritt mellom hver måling
w_d = delta_phi/T_s;  % vinkelhastighet i grader/s
noRounds = 10;         % antall runder rotasjon

% Registrerte vinkelmålinger i grader
VinkelPosMotorB = 0:delta_phi:360*noRounds;
noMeas = length(VinkelPosMotorB);   % antall målinger
t_k = 0:T_s:(noMeas-1)*T_s;   % tidsvektor

% Beregner sinusverdien av vinkelposisjonen
% gitt som y(t) = sin(w*t)

y = sind(w_d*t_k);    % sind() for degrees

% Plotter hvordan motoren har rotert
figure
set(gcf,'position',[285 400 750 500])
subplot(2,2,1)
plot(t_k,VinkelPosMotorB,'bo:')
xlabel('tid [s]')
ylabel('[$^{\circ}$]')
title('Vinkelposisjon $\varphi^{\circ}$ ved',...
    ['vinkelhastighet $\omega_d = ',num2str(w_d),...
    '~[^{\circ}$/s]'])

subplot(2,2,2)
plot(t_k,y,'b:o')
title('Estimat p{\aa} rotasjon, $y_k=\sin(\omega_d{\cdot} t_k)$')
xlabel('tid [s]')

uiwait(msgbox('Trykk OK for å oppsample'))

% Rekonstruerer signalet ved å legge inn
% n målinger og tidspunkt mellom hver
% samplet måling og måletidspunkt. Bruker interp-kommando.
n = 7;
VinkelPosMotorB_up = interp(VinkelPosMotorB,n);
t_k_up = interp(t_k,n);

% Beregner sinus på ny, w_d er den samme
y_up = sind(w_d*t_k_up);

% Plotter på nyi samme figur
subplot(2,2,3)
plot(t_k,VinkelPosMotorB,'bo:')
hold on
plot(t_k_up, VinkelPosMotorB_up,'r.-')
ylabel('[$^{\circ}$]')
title('Oppsamplet vinkelposisjon $\varphi_{up}^{\circ}$')
xlabel('tid [s]')

subplot(2,2,4)
plot(t_k_up, y_up,'rx-')
hold on
plot(t_k,y,'bo:')
xlabel('tid [s]')
title('Reell rotasjon, $y_k=\sin(\omega_d{\cdot} t_{k,up})$')

%}