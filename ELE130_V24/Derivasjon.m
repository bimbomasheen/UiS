function Sekant = Derivasjon(FunctionValues, Timestep, options)
arguments
    FunctionValues (1,3) double
    Timestep (1,1) double
    options.metode (1,:) char = 'Bakover'
end

if strcmp(options.metode,'Bakover')
    % fyll inn

elseif strcmp(options.metode,'Forover')
    % fyll inn

elseif strcmp(options.metode,'Senter')
    % fyll inn

else
    errordlg('Feil metode spesifisert')
    return
end

end
