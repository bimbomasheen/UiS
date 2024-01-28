% Ny oppstartsmappe, kopier inn mappestien etter kommandoen cd
% cd 

% Setter noen default-verdier for figurer
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',16)
set(0,'defaultTextFontSize',16)
set(0,'defaultLineLineWidth', 1);

% Setter rekkefølge på fargene i plot
farger = [0 0 1;  % blå
          1 0 0;  % rødt
          0 1 0;  % grønt            
          0 0 0]; % svart
set(0, 'DefaultAxesColorOrder', farger);


