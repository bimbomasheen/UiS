 %% -----------------------------------------------------
% a)
%
% Innhold:
%  - Vektor med tilfeldige tall
%  - Bruk av diff-funksjonen 
%  - Plotte ulike elementer 
%  - Finne maksimum og mininimum 
%
% Eksempler på nyttige Matlab-funksjoner
%  - rng
%  - randn
%  - num2str
%  - plot
%  - max
%  - min
%  - diff
%  - subplot
%  - title med bruk av num2str 
%  - title som deles på 2 linjer
%--------------------------------------------------

clear; close all; clc

seed = 1
rng(seed)
a = randn(1,45)

figure
% First subplot
subplot(2,1,1)

% Plot a with blue color, solid line, and 'x' markers
plot(a, 'b-x')
grid
hold on

% Title
title(sprintf('Signal $a$, hvor $a$(3) = %.2f', a(3)))

% Second subplot
subplot(2,1,2)

% Diff function on a
a1 = diff(a)

% Plot a1 with red color, solid line, and 'x' markers
plot(a1, 'r-x')
grid
hold on

% X-label and title
xlabel('Indeks $k$')
title({'Signal $a1$', 'Elementvis differanse i $a$'})

% Max difference
[a1_max, x_max] = max(a1)
plot(x_max, a1_max, 'b*', 'MarkerSize', 10)

% Min difference
[a1_min, x_min] = min(a1)
plot(x_min, a1_min, 'g*', 'MarkerSize', 10)

% Legend for second subplot
legend('Signal $a1$', 'St{\o}rste positive differanse', 'St{\o}rste negative differanse', 'Location','eastoutside')

% Return to the first subplot
subplot(2,1,1)

% Line for max positive difference
x_pos = [x_max, x_max+1]
plot(x_pos, a(x_pos),'b','LineWidth',3)

% Line for max negative difference
x_neg = [x_min, x_min+1]
plot(x_neg, a(x_neg),'g','LineWidth',3)

% Legend for first subplot
legend('Signal $a$', 'St{\o}rste positive differanse', 'St{\o}rste negative differanse', 'Location','eastoutside')

%% 

%% -------------------------------------------------------

% b)
%
% Innhold:
%  - Matrise og vektor med tilfeldige tall
%  - Uttrekk av forskjellige elementer
%  - Bruk av find-funksjonen til å finne utvalgte verdier
%  - LaTeX-tips
%
% Eksempler på nyttige Matlab-funksjoner
%  - length, numel og size
%  - disp og sprintf
%  - for-løkker og if-setninger
%  - plot, grid, hold, xlabel
%  - ones
%  - find
%  - title med LaTeX-kode for ø og å, {\o} og {\aa}
%  - legend med LaTeX-kode for ø og å, {\o} og {\aa}
%  - legend med styring av plassering
%----------------------------------------------------------

clear; close all; clc

%--------------------------------------------------
% Spesifiserer styrte tilfeldige tall med seed
% Lager matrisen B_mat
%--------------------------------------------------
seed = 3;
rng(seed)
B_mat = randn(3,45);

% Forskjellige kommandoer for å finne ulike størrelser fra matrisen B_mat.
k = length(B_mat); 
tot_elements = numel(B_mat);   
[len,bre] = size(B_mat); 

% Alternativ bruk av size
len = size(B_mat, 1);  % size kan gi kun lengde
bre = size(B_mat, 2);  % size kan gi kun bredde

disp(['length gir største dimensjon en vei: ', k ,' elementer'])
disp(['numel gir totalt antall: ',num2str(tot_elements),' elementer'])
str = sprintf(['size gir dimensjon [',num2str(len),' x ',num2str(bre),']\n']);
disp(str)

%--------------------------------------------------
% Lager vektor b som rad 2, alle kolonner fra B_mat
%--------------------------------------------------
b = B_mat(2, :);

% Plot de tilfeldige tallene b med blå strek og 'x' som markør
% Lager først indeksvektor
x = 1:length(b);

figure(1)
set(1,'position',[600 400 700 400])
plot(x, b, 'b-x')

% Siden du skal plotte flere ting, benytter vi "hold on".
% Du kan også bruke bare "hold", men en ny "hold" slår av holdingen 
hold on

% Kan også bruke "grid on", for å låse grid på. Ny "grid" slår av.
grid
title('Signal $b$ med markering p{\aa} visse element')
xlabel('Indeks $k$')


%-----------------------------------------------------------
% Uttrekk og presentasjon/plotting av elementer fra vektor b
%-----------------------------------------------------------

% Første element av b
b1 = b(1);   
disp(['Første element   = ',num2str(b1)])


% Siste element av b
b2 = b(end);   
disp(['Siste element    = ',num2str(b2)])


% De 4 første elementene 
b3 = b(1:4);   
disp(['4 første element = ',num2str(b3)])
% Lager tilhørende indeksvektor for plotting
x3 = 1:4;      
% Plot disse med en relativt stor 'x' som markør i grønn farge
plot(x3,b3,'gx','MarkerSize',10)


% De 4 siste elementene
b4 = b(end-3:end);     
disp(['4 siste element  = ',num2str(b4)])
disp(' ') % alternativ til \n i sprintf
% Lager tilhørende indeksvektor for plotting, bruk indeks k
x4 = (k-3):k;  
% Plot disse med en relativt stor '*' som markør i grønn farge
plot(x4,b4, 'g*' ,'MarkerSize',10)


% max-funksjonen finner indeks x5 og største verdi av b som vi kaller b5
[b5 , x5] = max(b);
disp(['Ved indeks k = ',num2str(x5),...
    ' er største verdi b = ',num2str(b5)])
% Plot inn største verdi med en relativt stor 'o' som markør i rød farge
plot(x5,b5, 'ro' ,'MarkerSize',10)


% min-funksjonen finner indeks x6 og minste verdi av b, b6
[b6, x6] = min(b);
disp(['Ved indeks i = ',num2str(x6),...
    ' er minste verdi b = ',num2str(b6)])
% Plot inn minste verdi med en relativt stor '*' som markør i rød farge
plot(x6,b6, 'r*' ,'MarkerSize',10)


% find-funksjonen gir indeksene x7 til alle verdier mellom b_min og b_max
b_max = -0.5;
b_min = -1;
x7 = find(b >= b_min & b < b_max);  
b7 = b(x7);    % plukker ut tilhørende element
disp(['Ved indeks i = [',num2str(x7),']'])
disp(['er verdiene av b = [',num2str(b7,3),']'])
disp(' ')

% Plot verdiene mellom b_min og b_max rett plass
% Bruk en relativt stor 'o' som markør i svart farge
plot(x7,b7,'ko','MarkerSize',10)

% Indikerer max- og min-verdiene med en vannrett svart stiplet strek
plot(b_min*ones(1, k),'k--')
plot(b_max*ones(1, k),'k--')

% sjekk oppgaveteksten hvordan legendene skal se ut
legend('Signal $b$', ...
    'De 4 f{\o}rste', ...
    'De 4 siste', ...
    'Maksimumsverdi', ...
    'Minimumsverdi', ...
    ['$b<', num2str(b_max), '$ og $b>', num2str(b_min), '$'], ...
    'Location','northwest','numcolumns',2);


% Alternativ måte å finne verdiene mellom b_min og b_max
for i = 1:length(b)
    if b(i) >= b_min && b(i) < b_max
        disp(['Ved indeks i = ', num2str(i), ...
            ' er verdien av b = ', num2str(b(i), 3)])
    end
end


%--------------------------------------------------
% Plot de tilfeldige tallene b i ny figur, figur 2
%--------------------------------------------------
figure(2)
set(gcf,'position',[400 100 700 400])

% Plasserer figuren slik at den ikke kommer rett over 
% den første. Flytt figuren der du vil og bruk kommandoen
%   >> get(gcf,'position')
% og kopier inn verdien i vektoren under.
plot(x,b,'b-x')
hold on
grid
title('Signal $b$')
xlabel('Indeks $k$.')


% Plukker ut annenhvert element, fra og med første.
b8 = b(1:2:end);
% Lager indeksvektor for plotting
x8 = 1:2:length(b);
% Plot disse med en relativt stor '.' som markør i rød farge
plot(x8, b8, 'r.', 'MarkerSize', 12)

% Plukker ut tredjehver element, starter fra andre element.
b9 = b(2:3:end);
% Lager indeksvektor for plotting
x9 = 2:3:length(b);
% Plot disse med en relativt stor 'o' som markør i grønn farge
plot(x9, b9, 'go', 'MarkerSize', 12)


legend('Signal $b$',...
    'Annenhvert element, fra f{\o}rste element',...
    'Tredjehver element, fra andre element',...
    'location','northwest')


%% 




%% --------------------------------------------------
% c)
% 
% Innhold:
%  - Vektor med tilfeldige tall sortert
%  - Bruk av statistike operasjoner
%  - Bruk av avrundingsfunksjoner
%  - Overordnet tittel i subplot
%
% Eksempler på nyttige Matlab-funksjoner
%  - sort
%  - abs
%  - mean, std og median
%  - ones
%  - round, fix, floor, ceil
%  - sum og prod
%  - sgtitle
%--------------------------------------------------


clear; close all; clc

% Absoluttveriden av tilfeldige tall, sortert fra høyest til lavest verdi 
seed = 7;
rng(seed)
c = sort(abs(randn(1,45))+0.3,'descend');
k = length(c);

% Plot de til sorterte tallene i stigende rekkefølge
figure
plot(c,'r-x')
hold on
title('Signal $c$ med statistiske verdier')
xlabel('Indeks $k$')

% Beregner middelverdi c1, standardavvik c2 og medianverdi c3,
% og plotter inn disse som hhv. linjer og punkt
c1 = mean(c);
c2 = std(c);
c3 = median(c);
x3 = find(c == c3);   % finn indeksen x3 til medianverdien ut fra c3

plot(c1*ones(1,k),'k--')      % middelverdi
plot(c1 + c2*ones(1,k),'k:')  % standardavvik over middelverdi
plot(c1 - c2*ones(1,k),'k:')  % standardavvik under middelverdi 
plot(x3, c3,'k*','MarkerSize',20)

% Skriver ut verdiene med redusert antall desimaler,
% ved bruk av disp og sprintf
disp(['Middelverdi = ',num2str(c1, '%.2f')])
disp(['Standardavvik = ',num2str(c2, '%.2f')])
str = sprintf(['Median = ',num2str(c3, '%.2f'),'\n']);
disp(str)

% Beregner, finner og indikerer verdiene innenfor ett std-avvik.
% Bruk kombinasjon av c1 og c2 til i find(..) til å finne indeksene x4. 
x4 = find(c>=c1-c2 & c<=c1+c2);
c4 = c(x4);
plot(x4,c4,'bo')

% Beregner %-andelen elementer innenfor ett standardavvik
andel = length(x4)/length(c) * 100;
disp(['Andel elementer innenfor standardavviket: ',num2str(andel,3),'%'])
disp(' ')

legend('Signal $c$',...
    ['Middelverdi  $\bar{c}=$',num2str(c1,2)],...
    ['Standardavviket over middelverdi, $\sigma$= ',num2str(c2,2)],...
     'Standardavviket under middelverdi',...
    ['Median-verdien = ',num2str(c3,2)],...
    'Elementene innenfor ett standardavvik',...
    'location','northeast')

%--------------------------------------------------
% Ulike avrundinger av c
%--------------------------------------------------
c5 = round(c);    % round
c6 = fix(c);    % fix
c7 = floor(c);    % floor
c8 = ceil(c);    % ceil

figure
set(gcf,'position',[400 100 700 400])
subplot(2,2,1)
plot(c,'r-x')
grid
hold on
plot(c5,'b-x')
title('round')

subplot(2,2,2)
plot(c,'r-x')
grid
hold on
plot(c6,'b-x')
title('fix')

subplot(2,2,3)
plot(c,'r-x')
grid
hold on
plot(c7,'b-x')
title('floor')

subplot(2,2,4)
plot(c,'r-x')
grid
hold on
plot(c8,'b-x')
title('ceil')

sgtitle('Forskjellige avrundinger av signal $c$')

% Beregn summen c7 og produktet c8 av alle elementene i c
c7 = sum(c);    % summen
disp(['Sum av alle element  = ',num2str(c7)])

c8 = prod(c);    % produktet
disp(['Produkt av alle element  = ',num2str(c8)])
%% 










%% --------------------------------------------------
% d)
% 
% Innhold:
%  - Bruk av for-løkker til rekursiv beregning.
%  - Sammenligning av for-løkker med elementvise operasjoner
%
% Eksempler på nyttige Matlab-funksjoner
%  - sqrt
%  - ^
%  - .^ (elementvise operasjoner)
%  - for-løkker
%  - set
%  - stem
%  - bar
%--------------------------------------------------

clear; close all; clc

% Vektoren x skal se være: x = 1.5, 2, 2.5, ..., 21
x = ..;
M = ..;  % antall elementer

%------------------------------------------
% Beregner d1 på 2 alternative måter
%------------------------------------------

% Alternativ 1
% Initialisering utenfor for-løkka
d1(1) = sqrt(..);
for k = ..:..
    d1(k) = ..;
end
disp(['d1(end) = ',num2str(d1(end),4)])
clear d1

% Alternativ 2
for k = ..:..
    % Initialisering innenfor for-løkka
    if k.. 
        d1(1) = sqrt(..);
    else
        d1(k) = ..;
    end
end

disp(['d1(end) = ',num2str(d1(end),4)])


%------------------------------------------
% Beregner d2 på en linje
%------------------------------------------
d2 = sum..sqrt..;
disp(['d2 = ',num2str(d2)])


%------------------------------------------
% Beregner d3
%------------------------------------------
d3(1) = ...;
for k=2:..
    d3(k) = ..;
end
disp(['d3(end) = ',num2str(d3(end))])

%------------------------------------------
% Beregner d4
%------------------------------------------
d4 = sum(x(..)..);
disp(['d4 = ',num2str(d4)])

% Plotter
figure
set(gcf,'position',[260 200 500 700])
subplot(3,1,1)
stem(x,'bx')
grid
title('$x(k)$')
xlabel('Indeks $k$')

subplot(3,1,2)
plot(d1,'rx')
grid
% lag tittelen riktig, sjekk oppgaven
title($d .. (k)  sum _ ^  sqrt  $')
xlabel('Indeks $k$')

subplot(3,1,3)
bredde = 0.4;
bar(d3,bredde,'g')
grid
% lag tittelen riktig, sjekk oppgaven
title($d ..  sum $')
xlabel('Indeks $k$')














%% 
% --------------------------------------------------- 
% e)                        
%
% Innhold:
%  - Analyse av rekkeutviklingen bak eksponentialfunksjonen
%  - Sammenligning mellom Matlab-funksjon og egenprodusert kode
%  - Sammenligning av resultat fra for-løkke og elementvise operasjoner
% 
% Eksempler på nyttige Matlab-funksjoner
%  - exp
%  - factorial
%  - log10 og log
%  - ./  (elementvise operasjoner)
%  - .*  (elementvise operasjoner)
%  - format
%  - sgtitle
%  - axis
%--------------------------------------------------

clear; close all; clc

x=2; 
M=40;

%------------------------------------
% Beregner e_xM1 som ett tall som øker
% for hver runde i for-løkka
%------------------------------------
e_xM1=0;
for n=..:..
    e_xM1 = ..
end
disp(['for-løkke:      e_xM1 = ',num2str(e_xM1)])

%------------------------------------
% Beregner e_xM2
%------------------------------------
e_xM2 = [];
n = ..:..;
e_xM2 = sum(..factorial ..);
disp(['elementvis:     e_xM2 = ',num2str(e_xM2)])

disp(['Avvik mellom for-løkke og elementvis = ',num2str(e_xM2 - e_xM1)])

%---------------------------------------------------
% Beregner e^x (kalles e_x) ved bruk av exp-funksjonen
%---------------------------------------------------
e_x = exp(..);
disp(['exp-funksjonen: e^x = ',num2str(e_x)])

disp(['Avvik mellom exp-funksjon og for-løkke = ',num2str(e_x-e_xM2)])

%------------------------------------
% Beregner avvik som funksjon av M
%------------------------------------
clear e_xM2  % sletter variabelen e_xM2
M = 0:100;
for k = 1:length(M)
    n = 0:..;
    e_xM2(k) = sum(.. factorial..);
    avvik(k) = e_x - e_xM2(k);
end

%------------------------------------
% Plotter avviket med lineær skala
%------------------------------------
figure(1)
set(gcf,'position',[400 600 1000 500])
subplot(1,3,1)
plot(M, avvik, 'b')
xlabel('{\O}vre grense $M$')
ylabel('avvik, line{\ae}r skala')
grid
axis([0 100 -1 10])

%------------------------------------
% Plotter avviket med logaritmisk skala
%------------------------------------
subplot(1,3,2)
semilogy(M, avvik, 'b')
xlabel('{\O}vre grense $M$')
ylabel('avvik, logaritmisk skala')
grid
axis([0 100 1e-16 10])

%------------------------------------
% Forskjell på log(e^x) og log10(e^x) for x=2
%------------------------------------
disp(' ')
disp(['Naturlige logaritmen log(e^x) = ',num2str(log(e_x))])
disp(['  10''er logaritmen log10(e^x) = ',num2str(log10(e_x))])

%------------------------------------
% Plotter 10'er-logaritmen av avviket
%------------------------------------
subplot(1,3,3)
avvik_log10 = log10(avvik); 
plot(M, avvik_log10, 'b')
xlabel('{\O}vre grense $M$')
ylabel('$\log_{10}$(avvik), verdi tilsvarer eksponenten')
grid
axis([0 100 -16 1])

% tittel over alle titler
sgtitle(['Avvik mellom $e^x$ og ',...
    '$e(x,M)$ for $x=2$ og {\o}kende $M$'])


%------------------------------------
% Presenterer avvikets størrelsen
%------------------------------------
disp(' ')
disp(['Siste element i log10(avvik) med ' ...
    '''format short'' og ''format long'':'])
verdi = avvik_log10(..);   % plukk ut riktig element
format short
verdi
format long
verdi
format short


%---------------------------------------------
% Presenterer grafisk hvert ledd i ligningen
% x^n * (1/n!) = ledd_1 * ledd_2
%---------------------------------------------
M = 15;
n = 0:M;
ledd_1 = ..;   % første ledd i ligningen
ledd_2 = ..;   % andre ledd i ligningen
uttrykk = ledd_1..ledd_2;  % elementvis produkt

for k=1:numel(n)
    % summen frem til k
    e_xM3(k) = sum(uttrykk(..:..));
end

figure
set(gcf,'position',[800, 350, 700, 550])
subplot(2,2,1)
plot(n, ledd_1, 'b-x')
xlabel('$n$')
legend('$x^n$ for $x{=}2$')
grid

subplot(2,2,2)
plot(n, ledd_2, 'r-x')
xlabel('$n$')
legend('$\frac{1}{n}!$')
grid

subplot(2,2,3)
plot(n, uttrykk ,'g-x')
xlabel('$n$')
legend('$x^n {\cdot}\frac{1}{n!}$ for $x{=}2$')
grid

subplot(2,2,4)
plot(n,e_xM3,'k-x')
xlabel('{\O}vre grense $M$')
legend('$e(x,M){=}\sum_{n=0}^M \frac{x^n}{n!}$ for $x{=}2$',...
    'location','best')
grid


%---------------------------------------------
% Presenterer leddene logaritmisk slik at det 
% er mulig å observere at: 
%    log10(a*b) = log10(a) + log10(b)
%---------------------------------------------
figure
set(gcf,'position',[370 320 700 550])
subplot(2,2,1)
plot(n,log10(ledd_1),'b-x')
grid
xlabel('$n$')
legend('$\log_{10}(x^n)$ for $x{=}2$',...
    'location','northwest')

subplot(2,2,2)
plot(n,log10(ledd_2),'r-x')
grid
xlabel('$n$')
legend('$\log_{10}(\frac{1}{n!})$',...
    'location','northeast')

subplot(2,2,3)
plot(n,log10(uttrykk),'g-x')
legend(['$\log_{10}(x^n {\cdot}\frac{1}{n!})$',...
    ' for $x{=}2$'],...
    'location','northeast')
grid
xlabel('$n$')

subplot(2,2,4)
plot(n,log10(ledd_1),'b-x')
hold on
grid
plot(n,log10(ledd_2),'r-x')
plot(n,log10(uttrykk),'g-x')
legend('$\log_{10}(x^n)$ for $x{=}2$',...
    '$\log_{10}(\frac{1}{n!})$',...
    'bl{\aa} + r{\o}d',...
    'location','best')
axis([0 15 -15 5])
xlabel('$n$')



%---------------------------------------------
% Presenterer leddene logaritmisk slik at det 
% er mulig å observere at: 
%    log10(a/b) = log10(a) - log10(b)
%---------------------------------------------

% Kopier koden for figuren rett over og gjør
% de små nødvendige justeringene i koden


%  ---------------------------------------------------------









%% 
% --------------------------------------------------- 
% f)                        
%
% Innhold:
%  - Bruk av histogram, styring av intervall
%  - Styring av akseverdier
%  - Plassere tekst i figurer
%
% Eksempler på nyttige Matlab-funksjoner
%  - histogram
%  - ceil, floor
%  - axis
%  - xlim
%  - text
%  - sprintf og ylabel
% --------------------------------------------------- 

clear; close all; clc

seed = 1;
rng(seed)        % styring av tilfeldige tall

% Signal x
n = 100;       % antall elementer i x-vektoren
x = abs(randn(1,n)*10 + 30);

%--------------------------------------------
% Plotter tallene i x-vektor i subplot(2,1,1)
%--------------------------------------------
subplot(2,1,1)
plot(x,'b:o')
ylabel('Tallverdi')
xlabel('indeks')
title('Tilfeldige tall $x$')

%---------------------------------------------------
% Plotter histogrammet til x i subplot(2,1,2)
% samtidig som vi henter ut x_prop som gir 
% mye informasjon om histogrammet, se Command Window.
% 
% I koden nedenfor kan du velge å benytte 
%    - x_prop.BinEdges 
%    - x_prop.BinWidth 
%    - x_prop.BinLimits 
% til å hente ut data fra x og plotte i subplot(2,1,1). 
% Alternativt kan du hardkode grenseverdiene.
%----------------------------------------------------
figure(1)
set(1,'position',[400 400 700 400])
subplot(2,1,2)
x_prop = histogram(x)  % ---> sjekk Command Window

% Tekststreng for ylabel. 
% Bruker sprintf for å kunne fordele teksten opp på 2 linjer
str = sprintf('Antall elementer i $x$ \n innenfor intervallene');
ylabel(str)
xlabel('Intervall i $x$')
title('Automatisk inndeling i intervall')

%------------------------------------------------
% Fremheving av elementene fra x-vektor i ett intervall
%   x_low < x < x_high
%      40 < x < 50
% Bruk enten xprop fra histogrammet 
% eller hardkodede verdier. 
%------------------------------------------------

% Velg riktig element BinEdges som gir x_min = 40
x_low = x_prop.BinEdges(5);
% x_low = 40;   % alternativt hardkodet

% x_high er ett intervall større
x_high = x_low + x_prop.BinWidth;

% Velger ut indeks-verdier fra intervallet
i = find(x>x_low  &  x<x_high);

subplot(2,1,1)
hold on 
plot(i, x(i), 'r*')
legend('Tilfeldige tall $x$',...
    ['Elementene i intervallet $',num2str(x_low),...
    '{<}x{<}$',num2str(x_high)],...
    'location','south')


%--------------------------------------
% Ny figur.
% Plotter histogrammet til x hvor vi
% styrer antall intervall (nbins)
%--------------------------------------
figure(2)
set(2,'position',[400 100 700 400])
subplot(2,1,1)
% velg selv antall intervall
nbins = 40; 
histogram(x,nbins)
ylabel(str)
xlabel('Intervall i $x$')
title(['Selvvalgt oppdeling i {\tt nbins=',num2str(nbins),'} intervall'])
xlim(x_prop.BinLimits)


%-----------------------------------------
% Alternativ måte å styre histogrammet på er
% å styrer intervallgrensene og 
% intervallbredden (edges) som
% fra:bredde:til
%------------------------------------------
subplot(2,1,2)

% Intervallgrensene, runder oppover og nedover
max_x = ceil(max(x));
min_x = floor(min(x));

% Velg intervallbredde
int_width = 2; 
edges = min_x:int_width:max_x;
histogram(x,edges)
ylabel(str)
xlabel('Intervall i $x$')
title(['Selvvalgt oppdeling av intervall som {\tt edges=[',...
    num2str(min_x),':',num2str(int_width),':',num2str(max_x),']}'])
xlim(x_prop.BinLimits)
hold on

%-----------------------------------------------------------
% Skissering av middelverdi og standardavvik i histogrammet
%-----------------------------------------------------------
mean_value = mean(x);
std_value = std(x);

% x-koordinaten til text-funksjonen
x1 = mean_value;

% Benytter y-akseinformasjon for å både plotte loddrett 
% strek samt å plassere tekst 
ax = axis;    % leser av aksene
y1 = ax(3);  % plukker ut laveste y-akseverdi
y2 = ax(4);  % plukker ut høyeste y-akseverdi


% Definerer fargen darkgreen ut fra [R G B], 
% hvor [1 1 1] = hvitt og [0 0 0] = svart
darkgreen = [0 0.5 0.0];  

% Plotter loddrett strek i farge darkgreen
plot([x1 x1],[y1 y2],'color',darkgreen, 'LineWidth',3)

% Plasserer inn en tekststreng hvor vi hardkoder justering på 
% plasseringen med å gange med en tallverdi større eller mindre enn 1
text(x1*1.1, y2*0.9,...
    ['Middelverdi $\bar{x}$=',num2str(mean_value,3)],...
    'color',darkgreen)

% x- og y-koordinatene til text
x1 = mean_value;  % gjentar denne
x2 = mean_value + std_value;
y1 = (ax(4)+ax(3))/1.5;

% Definerer fargen darkred ut fra [R G B]
darkred = [0.9 0.2 0.2];

% Plotter vannrett strek i farge darkred
plot([x1 x2],[y1 y1],'color', darkred,'LineWidth',3)
text(x1*1.05, y1*1.05 ,...
    ['Standardavvik $\sigma$=',num2str(std_value,3)],...
    'color',darkred)
% ---------------------------------------------------------













%% 
% --------------------------------------------------- 
% g)   
%
% Innhold:
%  - Bruk av findpeaks-funksjonen
%
% Eksempler på nyttige Matlab-funksjoner
%  - sin
%  - findpeaks
%  - text
%  - arrow (hentet fra Mathworks File Exchange)
%    https://se.mathworks.com/matlabcentral/fileexchange/278-arrow
% --------------------------------------------------- 


clear; close all; clc

delta_t = 0.02;
t = 0:delta_t:5;
u = 3.2*sin(4*t)+2;

% Det er to måter å bruke findpeaks på:

% --------------------------------------------------- 
% 1) Uten returargument plottes figur med toppunkter, 
%    enten som funksjon av indeksene eller 
%    som funksjon av en definert tidsvektor t.
% --------------------------------------------------- 
figure
subplot(2,1,1)
findpeaks(u)  
xlabel('indeks $k$') 
title('Kommandoen {\tt findpeaks(u)} plotter $u(k)$')

subplot(2,1,2)
findpeaks(u, t)
xlabel('tid [s]') 
title('Kommanoden {\tt findpeaks(u,t)} plotter $u(t)$')

% --------------------------------------------------- 
% 2) Med bruk av returargument plottes ikke figuren,
%    men toppunkt-verdiene og indeksene returneres,
%    her kalt u_pks og u_locs
% --------------------------------------------------- 
[u_pks,u_locs] = findpeaks(u);

% Tidspunkt for maksverdi 
t_pks = t(u_locs);

% skriver ut til Command Window
disp(['Ved t = [',num2str(t_pks,3),']'])
disp([' er u = [',num2str(u_pks,3),']'])


%-------------------------------------------------------
% Automatisk beregning av periode, vinkelfrekvens,
% likevektsverdier og amplituder.
% Plotter først signalene u(t) og y(t) med findpeaks-funksjonen
%-------------------------------------------------------
y = 1.4*sin(4*t-1.9)+3;

% Plotter toppene i u(t) og y(t)
figure
findpeaks(u,t)
hold on
findpeaks(y,t)


% Finner maksverdiene og tilhærende indeks
[u_pks,u_locs] = findpeaks(u);
[y_pks,y_locs] = findpeaks(y);

% Finner tidspunktet tilsvarende indeksene
t_pks_u = t(u_locs);
t_pks_y = t(y_locs);

% Beregner T_p og omega
T_p = mean(diff(t_pks_y));
w = 2 * pi / T_p;

% Beregner likevekstverdi/konstantleddet
u_A = mean(u);
y_A = mean(y);

% Beregner amplitudene
U = (max(u)-min(u)) / 2;
Y = (max(y)-min(y)) / 2;

% Beregner faseforsyvningen
delta_t_phi = t_pks_u(1) - t_pks_y(1);
phi = w * delta_t_phi;
phi_d = phi * (180 / pi);


%-------------------------------------------------------
% Legger på piler og linjer og tekst som viser
% beregnede verdier.
%-------------------------------------------------------

% Periode T_p og vinkelfrekvens w
xy1 = [t_pks_u(1),u_pks(1)];
xy2 = [t_pks_u(2),u_pks(2)];
arrow(xy1,xy2,'Ends',[1,2])
text(0.8,5.7,['$T_p{\approx}$',num2str(T_p,3),' s, $\Rightarrow ',...
' \omega{\approx}',num2str(w,3),'$ rad/s'])

% Tidsforskyvning delta_t_phi og faseforskyvning phi
xy1 = [t_pks_u(1),1];
xy2 = [t_pks_y(1),1];
arrow(xy1,xy2,'Ends',[1,2])
text(0.2,0.6,['$\Delta t_{\varphi}{\approx}',num2str(delta_t_phi),'$ s'],...
    'backgroundcolor','white')
text(0.5,0.1,'$\Downarrow$','backgroundcolor','white')
text(0.2,-0.4,['$\varphi{\approx}',num2str(phi,3),'$'],...
    'backgroundcolor','white')
text(0.5,-0.9,'$\Downarrow$','backgroundcolor','white')
text(0.2,-1.3,['$\varphi^{\circ}{\approx}',num2str(phi_d,4),'^{\circ}$'],...
    'backgroundcolor','white')

%   Første loddrette strek
x = [t_pks_u(1),t_pks_u(1)];
y = [1,u_pks(1)];
plot(x,y,'k:')
%   Andre loddrette strek
x = [t_pks_y(1),t_pks_y(1)];
y = [1,y_pks(1)];
plot(x,y,'k:')

% Likevektslinjen u_A
x = [0,t(end)];
y = [u_A,u_A];
plot(x,y,'b--')
text(2.5,u_A,['$u_A{\approx}$',num2str(u_A)],'backgroundcolor','white')

% Likevektslinjen y_A
x = [0,t(end)];
y = [y_A,y_A];
plot(x,y,'r--')
text(1.2,y_A,['$y_A{\approx}$',num2str(y_A)],'backgroundcolor','white')

% Amplituden U 
xy1 = [t_pks_u(3),u_A];
xy2 = [t_pks_u(3),u_pks(3)];
arrow(xy1,xy2,'Ends',1,'width',1)
text(3.3,3.5,['$U{\approx}$',num2str(U)],'backgroundcolor','white')

% Amplituden Y 
xy1 = [t_pks_y(2),y_A];
xy2 = [t_pks_y(2),y_pks(2)];
arrow(xy1,xy2,'Ends',1,'width',1)
text(2.2,3.5,['$Y{\approx}$',num2str(Y)],'backgroundcolor','white')

legend('$u(t)=3.2\sin(4t)+2$',...
    'max-verdier i $u(t)$',...
    '$y(t)=1.4\sin(4t-1.9)+3$',...
    'max-verdier i $y(t)$',...
    'location','southeast')

title('Signalene $u(t)$ og $y(t)$ med automatisk beregnede verdier')














%% --------------------------------------------------------- 
% h)   
%
% Innhold:
%  - Bruk av pause til å demonstrere utvikling i for-løkker
% 
% Eksempler på nyttige Matlab-funksjoner
%  - sin
%  - pause
%  - ctrl C
%  - hold on/hold off
%  - tic og toc
% --------------------------------------------------------- 

clear; close all; clc

y1 = 0;
y2 = 0;
t=0:0.01:20;
w=2;  % omega, rad/s

% Starter stoppeklokka, UTEN returvariabel. 
% For å få tak i medgått tid brukes toc direkte
tic 

% ------------------------------------------------
% Spesifiser M, lek gjerne litt med denne verdien
M = 4;    
% ------------------------------------------------

for n=1:M
    y1 = y1 + ..
    y2 = y2 + ..

    subplot(2,1,1)
    plot(t,y1,'b')
    % hold on     % test med og uten

    subplot(2,1,2)
    plot(t,y2,'r')
    % hold on     % test med og uten

    xlabel('tiden $t$')
    disp(['Runde ',num2str(n),' av ',num2str(M)])
    disp('Pause. Mulig du må trykke en tast.')

    % Starter ny stoppeklokke MED returvariabel "NyKlokke"
    % Stoppeklokken starter på ny for hver runde
    NyKlokke = tic;
    pause  % alternativt pause(..)

    % Avleser medgått tid for den nye stoppeklokken.
    % Argumentet til toc() er returvariabelen.    
    disp(['Pausen varte i ',num2str(toc(NyKlokke)),' sekund'])

    % Henter ut mellomtid med bruk av toc direkte. 
    % Tiden er fra stoppeklokken startet før for-løkken
    disp(['Medgått tid siden start er ',num2str(toc),' sekund'])
    disp('-----------------------------------------')
end

disp(' ')
disp(['Hele for-løkken brukte ',num2str(toc),' sekund'])
disp(' ')

subplot(2,1,1)
title(['$y_1(t) = \sum_{n=1}^{M}\frac{1}{2n-1}',...
    '\sin\bigl((2n-1){\cdot}\omega {\cdot}t\bigr)$',...
    ', $M{=}$',num2str(M),' $\omega{=}$',num2str(w)]);
subplot(2,1,2)
title(['$y_{2}(t) = \sum_{n=1}^{M} \frac{1}{n}',...
    '\sin(n{\cdot}\omega{\cdot} t)$',...
    ', $M{=}$',num2str(M),' $\omega{=}$',num2str(w)]);












%% --------------------------------------------------- 
% i)   
%
% Innhold:
%  - Bruk av area-funksjonen
%
% Eksempler på nyttige Matlab-funksjoner
%  - area
%  - cos
% --------------------------------------------------- 

clear; close all; clc

U = 4;
w = 7;
C = 0;

figure(1)
set(gcf,'position',[500 100 450 600])

% Juster på denne
ant_runder = 3;

for i = 1:ant_runder
    t = 0:0.01:5;
    u = U*sin(w*t)+C;
    y = (-U/w) * (cos(w * t) - 1) + C * t;

    subplot(3,1,1)
    plot(t,u,'Linewidth',2)
    hold on
    grid
    title('$u(t) = U*sin(w*t)+C$')

    subplot(3,1,2)
    area(t,u)
    grid
    hold on
    title('Fargelagt areal under $u(t)$')

    subplot(3,1,3)
    plot(t,y,'LineWidth',2)
    grid
    hold on
    title('$y(t) = (-U/w) * (cos(w * t) - 1) + C * t$')

    % justering av C og/eller vinkelfrekvens w
    % C = C - 0.4;
    w = w - 2;
end

t_avlest = 1.78
y_avlest = 0.71523











%% 
% --------------------------------------------------- 
% j)   
%
% Innhold:
%  - Fjerne uteliggere, konstantledd og deler av et signal
%  - Relevant for Lego-prosjektet
% 
% Eksempler på nyttige Matlab-funksjoner
%  - randi, rand og randn
%  - zeros
%  - NaN
%  - findpeaks med MinPeakDistance-argument
% --------------------------------------------------- 

clear; close all; clc

% Tilfeldige parametre til sinussignalet u(t) gitt som
%  -> u(t) = U sin(wt) + C 
% For å beregne w, U og C tas snittet av 30 rand()-tall
% som deretter justeres opp med en faktor

seed = 1;
rng(seed);
w = 3*mean(rand(1,30));   % tilfeldig uniformfordelt vinkelfrekvens rad/s
U = 5*mean(rand(1,30));   % tilfeldig uniformfordelt amplitude 
C = 3*mean(rand(1,30));   % tilfeldig uniformfordelt konstantverdi

% Rent sinussignal u_sin med tidsvektor t, 29 sekund lang
t = 0:0.1:29;
u_sin = U*sin(w*t);

% Lager forsinket sinussignal ved å sette inn
% k stykk 0'ere i begynnelsen av u_sin og samtidig  
% fjerner like mange elementer fra u_sin i slutten
% 
% For at variasjonen ikke skal være så stor ved 
% forskjellige rng-verdier, tas gjennomsnittet av 10
% randi()-tall for å beregne forsinkelsen k.
k = 6*round(mean(randi(10,1,10))); 
u_delay = [zeros(..) u_sin(..)];

% Legger på litt normalfordelt målestøy,
% hvor støyamplituden er lik 10% av amplituden U
v = 0.1*U*randn(1,length(t)); 

% Lager u(k) som forsinket sinus + konstant + målestøy
u = u_delay + C + v;   

% Lager to uteliggere
u(end) = max(u)*2;            % Uteligger 1
u(round(end/2)) = - max(u);   % Uteligger 2

figure
set(gcf,'position',[50 100 600 850])
subplot(4,1,1)
plot(t,u,'b')
grid
title('Signal $u(t)$ som funksjon av tiden $t$')
xlabel('tid [s]')

subplot(4,1,2)
plot(u,'b')  % kun ett argument i plot-kallet = plotting mot indeks
grid
title('Signal $u(k)$ som funksjon av indeks $k$')
xlabel('indeks $k$')

%--------------------------------------------------------
% Estimerer startindeks k_est for sinus i subplot(4,1,2)
%--------------------------------------------------------
k0 = ..; 

% Fjerner de første elementene i t og overskriver i ny t.
% Samme for u
t = t(..);
u = u(..);

% Justerer tidssignalet slik at t(1)=0
t = t - ..;

% Fjerner uteliggere
u_max = ..;
u_min = ..;
i = find(u>.. | u<..);
u(i) = NaN;

% Finner estimat av likevektslinjen
C_est = max(u).. min(u)..

% Justert signal som er den "rene" sinusbølgen
u = u - C_est;

subplot(4,1,3)
plot(t,u,'b')
grid
xlabel('tid [s]')

%---------------------------------------------------------------
% Prøve- og feilemetoden for å finne amplitude U og frekvens w.
% Modell 1.
%---------------------------------------------------------------

U_est1 = 1;
w_est1 = 1;
u_est1 = U_est1*sin(w_est1*t);

hold on
plot(t,u_est1,'r')
title('Justert signal $u(t)$, modell funnet manuelt')
legend('Justert signal $u(t)$',...
    ['Modell 1: $u(t)=',num2str(U_est1),...
    '\sin(',num2str(w_est1),'{\cdot}t)$'],...
    'Location','east')
axis('tight')


%----------------------------------------------------------
% Automatisk beregning av amplitude U og frekvens w. 
% Modell 2.
%----------------------------------------------------------
subplot(4,1,4)

% Her må du leke deg litt med forskjellige verdier for 
% MinPeakDistance og studere effekten av de
findpeaks(u,t,'MinPeakDistance',..)
[u_pks,t_locs] = findpeaks(u,t,'MinPeakDistance',..);


% Beregn amplitude U_est2 på alternative måter
%   - alternativ 1 ut fra største og minste verdi i u-signalet
%   - alternativ 2 som snitt av u_pks
U_est2 = (max .. min )/2;    % alternativ 1  
U_est2 = mean(u_pks);          % alternativ 2

% Vinkelfrekvens w_est
% Beregner først snittet av periodetid T_p 
T_p = mean(diff(t_locs)); 
w_est2 = ..;

% lager signalet u_est2
u_est2 = U_est2*sin(w_est2*t);

hold on
plot(t,u_est2,'r')
legend('Justert signal $u(t)$',...
    ['Toppunkter, $T_p \approx $',num2str(T_p),' s'],...
    ['Modell 2: $u(t)=',num2str(U_est2,3),...
    '\sin(',num2str(w_est2,4),'{\cdot}t)$'],...
    'Location','east')
title('Justert signal $u(t)$, modell funnet automatisk')
xlabel('tid [s]')
% ---------------------------------------------------------


