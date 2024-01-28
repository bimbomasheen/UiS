%% Celle nr 1 
clear; clc; close all

a = 1;
disp(a)

%% Celle nr 2, feil i denne gir problem
clear; clc; close all

% les i kompendiet, kap 2.2 side 14
b =  % ufullstending kode


%% Presentasjon av beregningsresultat
clear;clc

% vise en variabel 
a = 314.159265;
disp(a)   

% lage strengen str_a
str_a1 = ['a1 = ',num2str(a)];
disp(str_a1)

% bygge opp strengen inne i disp
disp(['a2 = ', num2str(a)])

% formatert tekst
str_a3 = sprintf('a3 = %f',a);
disp(str_a3)

% formatert tekst, styring av format
str_a4 = sprintf('a4 = %1.2f',a);
disp(str_a4)

% styring av format med num2str
str_a5 = ['a5 = ',num2str(a,4)];
disp(str_a5)

str_a6 = ['a6 = ',num2str(a,'%1.2f')];
disp(str_a6)

disp(['Element 4 er: ',str_a6(4)])

%% Presentasjon av beregningsresultat, del 2
clear;clc

a = 6;
a = a + 5;
a = a - 3;

% Modulo-divisjon
b = rem(a,3);  % rem for remainder
disp(['b = ',num2str(b)])  

% liggende vektor, enten komma eller ikke
c = [1 2 3 4];
disp(['c = ',num2str(c)]);

% transponert ved bruk av apostrof '
% d blir staaande vektor
d = c';  
% disp(['d = ',num2str(d)]);  % gir feil
disp('d = ')
disp(d)


%% Indeksering, første indeks og initialverdi
clear;clc

x(1) = 3;
for k = 2:5
   x(k) = x(k-1) + 1;
end
disp(x)


%% Indentering, vil fungere
clear;clc

  x(1) = 3;
for k = 2:5
x(k) = x(k-1) + 1;
  end
disp(x)


%% Datatyper og kommentarer
clear;clc

%{ 
Lange 
kommentarer 
kan 
lages slik
%}

% Datatyper
c = 3;        % c blir double
d = 6.4;      % d blir double
e = 5 + 3j;   % e blir double
f = int8(d);  % e blir int8
disp(['c = ',class(c)]) 
disp(['d = ',class(d)]) 
disp(['e = ',class(e)]) 
disp(['f = ',class(f)]) 

% Funksjonaliten med _a
% finnes ikke 
% i Matlab

% Matlab kan bruke både "j", "1j", "i" og "1i"
% som kompleks variabel. 
i1 = 5 + 3j;
i2 = 5 + 3*1j;
i3 = 5 + 3i;
i4 = 5 + 3*1i;
g = 5;
h = 3;
i5 = complex(g,h);

disp(i1);disp(i2);disp(i3);disp(i4);disp(i5)


%% For-løkker
clear;clc

% metode 1, dynamisk utvidelse av A
A = []; 
for k=1:10
    A(k) = k-1;
end
disp(['A = ',num2str(A),' ',class(A)])


% metode 2, preallokering med zeros
B = zeros(1,10); 
for k=1:10
    B(k) = k-1;
end
disp(['B = ',num2str(B),' ',class(B)])

% metode 3, preallokering med ones
% hvor for kort preallokering 
% gir INGEN feilmelding
C = ones(1,5); 
for k=1:10
    C(k) = k-1;
end
disp(['C = ',num2str(C),' ',class(C)])


%% Rekursive beregninger i for-løkker
clear;clc

u = 1:10;

% versjon 1)
y1 = sum(u);
disp(['y1 = ',num2str(y1)])


% versjon 2)
M = length(u);
y2 = u(1);  % initialisering
for k = 2:M
    y2 = y2 + u(k);
end
disp(['y2 = ',num2str(y2)])


% versjon 3)
y3(1) = u(1);  % initialisering
for k = 2:M
    y3(k) = y3(k-1) + u(k);
end
disp(['y3 = ',num2str(y3)])



%% If-setninger og logiske betingelser
clear; clc

a = input('Skriv inn tallet "a" mellom 0 og 10: ');

if a >= 5
    disp('  a >= 5 \n')    
elseif a > 0 && a < 5
    disp('  a > 0 og a < 5 \n')    
else
    disp('  a <= 0 \n')    
end

b = input('Skriv inn tallet "b" mellom 0 og 10: ');

if a == 1 || b == 9
    disp('  a == 1 eller b == 9 \n');
elseif a ~= 1 && b ~= 9
    disp('  a ~= 1 og b ~= 9 \n');
end



%% Plotting, tittel og legender
clear; close all; clc

Y1 = 2;
Y2 = 0.5;
t = linspace(0,2*pi,100);
y1 = Y1*sin(t);
y2 = Y2*cos(t);
grid
plot(t,y1,'b')
hold on
plot(t,y2,'r:')
plot(t,y1+y2,'g--')
title(['St{\o}rste tid = ',num2str(max(t),3)])
ylabel('y-verdi')
legend(['$y_1$, ampl. Y1 = ',num2str(Y1)],...
        ['$y_2$, ampl. Y2 = ',num2str(Y2)],...
        'y1+y2, bare tekst')
text(1,-1,'Tekst')
xlabel('tiden $t$ [s]') 


%% Deklarering av lister/array/vektorer/matriser
clear; clc

C1 = [1, 2, 3];  % eller C1 = [1 2 3]; 

% vektor = [start:step:slutt];
C2 = 4:2:8; 

% linspace(start,slutt,#verdier)
% Tredje parameter er valgfri. 
% Hvis utelatt så blir 
% 100 verdier lagt inn
C3 = linspace(2,8,4);

% Initialiser en 2D array/matrise.
C4 = [1 2 3; 4 5 6];

disp(['C1 = ',num2str(C1)])
disp(['C2 = ',num2str(C2)])
disp(['C3 = ',num2str(C3)])
disp('C4 = ')
disp(C4)

% Skriv ut elementet med tallet 5
% fra matrisen C4
disp(['C4(2,2) = ',num2str(C4(2,2))])

% Skriv ut rekke 1 og 2 i 
% kolonne 1 av matrisen C4.
% Skriver ut tallene 1 og 4
for i = 1:2
    str = sprintf('C4(%d,1) = %d',i,C4(i,1));
    disp(str)
end
disp(' ')



%% Uttrekk av enkeltelement fra en liste/vektor
clear;clc

% 2 alternative måter
x = 1:10;
x = linspace(1,10,10);
disp(['x = ',num2str(x)])

disp('Første element')
disp(['    x(1) = ',num2str(x(1))])

disp('Nest siste element')
disp(['x(end-1) = ',num2str(x(end-1))])
disp(['    x(9) = ',num2str(x(9))])

disp('Siste element')
disp(['  x(end) = ',num2str(x(end))])
disp(['   x(10) = ',num2str(x(10))])

disp('Dersom x stadig blir større,') 
disp('bruker k som argument')
k = length(x);
disp(['k = ',num2str(k)])

disp('Nest siste element')
disp(['  x(k-1) = ',num2str(x(k-1))])

disp('Siste element')
disp(['    x(k) = ',num2str(x(k))])



%% Slicing, uttrekk av flere element fra en liste/vektor
clear;clc

x = 1:10;
disp(['x = ',num2str(x)])

disp('Første tre element')
disp(['        x(1:3) = ',num2str(x(1:3))]); 

disp('Andre og tredje element')
disp(['        x(2:3) = ',num2str(x(2:3))]); 

disp('Tredje element')
disp(['          x(3) = ',num2str(x(3))]);

disp('Tredje siste og nest siste')
disp(['x(end-2:end-1) = ',num2str(x(end-2:end-1))])
disp(['        x(8:9) = ',num2str(x(8:9))])

disp('De to siste elementene')
disp(['  x(end-1:end) = ',num2str(x(end-1:end))])
disp(['       x[9:10] = ',num2str(x(9:10))])

k = 6;
disp(['Spesifiserer vilkårlig k = ',num2str(k)])

disp('Plukker ut elementer')
disp(['    x(k-1:end) = ',num2str(x(k-1:end))])
disp(['      x(k-1:k) = ',num2str(x(k-1:k))])


disp('Indeksere utover lengden til x')
disp('gir feil ved slicing')
k = length(x);
disp(['Antall elementer i x, k = ',num2str(k)])
disp(['  x(k-1:k+1) = ',num2str(x(k-1:k+1))])




%% Stepping i en liste/vektor
clear;clc

x = 1:10;
disp(['x = ',num2str(x)])
% x[start:step:stop]

disp('Fra første element,') 
disp('til og med element Y,') 
disp('for hvert andre element.')
disp(['    x(1:2:end) = ',num2str(x(1:2:end))])

disp(['  x(1:2:end-1) = ',num2str(x(1:2:end-1))])

disp(['  x(1:2:end-2) = ',num2str(x(1:2:end-2))])

disp('Fra fjerde siste, til siste,')
disp('for hvert andre element')
disp(['  x(end-3:2:end) = ',num2str(x(end-3:2:end))])

disp('Fra indeks 6 til indeks 10,')
disp('for hvert andre element')
disp(['       x(6:2:10) = ',num2str(x(6:2:10))])


disp('Dersom x stadig blir storre, ')
disp('bruker k som argument.')
k = length(x);
disp(['Antall elementer i x, k = ',num2str(k)])

disp('Tredje siste og siste element')
disp(['    x(k-2:2:k) = ',num2str(x(k-2:2:k))])



disp('Tilsvarende funksjonalitet')
disp('finnes ikke i MATLAB')
disp('og gir feilmelding ')


%% Vilkårlige element i en liste
clear;clc

x = 1:10;
disp(['x = ',num2str(x)])

element = [4,8,10];
x1 = x(element);
disp(['   x([4,8,10]) = ',num2str(x1)])

% tilsvarende fumksjonalitet finnes 
% ikke i Matlab
disp(' ')

% dersom x stadig blir storre, 
% bruker k som argument
k = length(x);

element = [k-6,k-2,k];
x_k = x(element);
disp(['x([k-6,k-2,k]) = ',num2str(x_k)])



%% Vilkårlige element og vektorer fra en matrise 
clear; clc
x = 1:12;
A = reshape(x,3,4); % lager matrise av vektoren x
A = A'; % transponserer (roterer den om diagonalen) matrisen 
A
rad2 = A(2,:)     % rad 2, alle kolonner
kol3 = A(:,3)     % alle rader, kolonne 3
element = A(3,3)  % rad 3, kol 3
vektor = A(3,1:2) % rad 3, kolonne 1 til 2


%% Egne funksjoner
clear;clc

x = 2;
y = 3;
z = myfunc(x,y);
disp(['z = ',num2str(z)])  % skriver ut tallet 8 = 2^3 

% funksjonen ligger helt nederst i skriptet

%% Egne funksjoner
clear;clc

x = 2;
y = 3;
[z1, z2] = myfunc2(x,y);
disp(['z1 = ',num2str(z1)])  % skriver ut tallet 8 = 2^3 
disp(['z2 = ',num2str(z2)])  % skriver ut tallet 8 = 2^3 

% funksjonen ligger helt nederst i skriptet


%% Normalfordelte tilfeldige tall
clear;close all;clc

% styring av tilfeldighetene
seed = 15;
rng(seed);

% Normalfordeling
x1 = randn(1,1000);
x1_mean = mean(x1);
x1_std = std(x1);

figure
set(gcf,'Position',[560 400 550 430])
% Plotting av x 
subplot(2,1,1)
plot(x1)
title(['Tilfeldig signal $x_1$ med {\tt randn}-funksjonen,',...
    ' seed = ',num2str(seed)])

subplot(2,1,2)
bins = 50;
histogram(x1,bins)
xlabel(['Tallverdiene i $x_1$ delt opp ',num2str(bins),' intervall'])
ylabel(['Antall elementer'])
title('Fordelingen med styrt oppdeling')



%% Uniformfordelte tilfeldige tall
clear;close all;clc

% styring av tilfeldighetene
seed = 15;
rng(seed);

% Uniformfordeling
x2 = rand(1,1000);
x2_mean = mean(x2);
x2_std = std(x2);

figure
set(gcf,'Position',[560 400 550 430])
% Plotting av x 
subplot(2,1,1)
plot(x2)
title(['Tilfeldig signal $x_2$ med {\tt rand}-funksjonen,',...
    ' seed = ',num2str(seed)])

subplot(2,1,2)
bins = 50;
histogram(x2,bins)
xlabel(['Tallverdiene i $x_2$ delt opp ',num2str(bins),' intervall'])
ylabel(['Antall elementer'])
title('Fordelingen med styrt oppdeling')



%% Uniformfordelte tilfeldige heltall

clear;close all;clc

% styring av tilfeldighetene
seed = 15;
rng(seed);

% Uniformfordeling av heltall
max_verdi = 5;
x3 = randi(max_verdi,1,100);
x3_mean = mean(x3);
x3_std = std(x3);

figure
set(gcf,'Position',[560 400 550 430])
% Plotting av x 
subplot(2,1,1)
plot(x3,'x-')
axis([0 100 0 6])
title(['Tilfeldig signal $x_3$ med {\tt randi}-funksjonen,',...
    ' seed = ',num2str(seed)])

subplot(2,1,2)
bins = 50;
histogram(x3,bins)
xlabel(['Tallverdiene i $x_3$ delt opp ',num2str(bins),' intervall'])
ylabel(['Antall elementer'])
title('Fordelingen med styrt oppdeling')




%% Litt om statistiske funksjoner
clear;close all;clc

% illustrasjon på varians
x = [6.1, 6.4, 4.1, 2.7, 6.3, 1.2, 5.3, 1.1, 0.9, 2.2];

x_mean = mean(x);
sigma2 = var(x);
sigma = std(x);

figure
set(gcf,'Position',[560 200 700 750])

subplot(3,1,1)
plot(x,'b*',MarkerSize=10)
xlabel('indeks $k$')
hold on
plot(ones(1,length(x))*x_mean,'r*:',MarkerSize=10)
axis([0 length(x) 0 7])
for i=1:length(x)
    plot([i i],[x(i),x_mean],'k-',LineWidth=1)
    text(i+0.1,x_mean+(x(i)-x_mean)/2,...
        [num2str((x(i)-x_mean)^2,2)],...
        'interpreter','latex')
end
legend('datasett $x$',...
    'middelverdi $\bar{x}$',...
    '$(x(k)-\bar{x})^2$',...
    'interpreter','latex')
title('Datasettet $x$ med indikering av tallverdiene $(x(k)-\bar{x})^2$')


subplot(3,1,2)
plot(x,'b*',MarkerSize=10)
xlabel('indeks $k$')
hold on
x_mean = mean(x);
plot(ones(1,length(x))*x_mean,'r*:',MarkerSize=10)
axis([0 length(x) 0 7])
abs_std = 0;
for i=1:length(x)
    plot([i i],[x(i),x_mean],'k-',LineWidth=1)
    text(i+0.1,x_mean+(x(i)-x_mean)/2,...
        [num2str(abs(x(i)-x_mean),2)],...
        'interpreter','latex')
    abs_std = abs_std + abs(x(i)-x_mean);
end
abs_std = abs_std/length(x);
legend('datasett $x$',...
    'middelverdi $\bar{x}$',...
    '$|x(k)-\bar{x}|$',...
    'interpreter','latex')
title('Datasettet $x$ med indikering av tallverdiene $|x(k)-\bar{x}|$')


subplot(3,1,3)
histogram(x,9)
title(['Fordelingen av $x$ med $\bar{x}=$',num2str(x_mean),...
    ', $\sigma^2 = $',num2str(sigma2,3),...
    ', $\sigma = $',num2str(sigma,3)])
xlabel('Intervall for verdier i $x$')
axis([0 7 0 4])
hold on
plot([x_mean x_mean], [0 4],'r')
plot([x_mean-sigma x_mean], [2 2],'b')
plot([x_mean x_mean+sigma], [1.5 1.5],'b')
legend('fordeling','middelverdi',...
    ['bredde tilsvarende $\sigma = $',num2str(sigma,3)],...
        'interpreter','latex')


%% Egne funksjoner

% funksjonen må ligge nederst i .m-filen
function Value = myfunc(arg1,arg2)
    % ^ er eksponential-operator
    % Bruker lokale variable
    % MÅ definere egen returvariabel
    Value = arg1^arg2; 
end

function [Val1, Val2] = myfunc2(arg1,arg2)
    % ^ er eksponential-operator
    % Bruker lokale variable
    % MÅ definere egen returvariabel
    Val1 = arg1^arg2; 
    Val2 = arg1+arg2; 
end
