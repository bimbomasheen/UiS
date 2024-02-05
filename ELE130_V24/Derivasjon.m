function Sekant = Derivasjon(FunctionValues, Timestep, options)
arguments
    FunctionValues (1,3) double
    Timestep (1,1) double
    options.metode (1,:) char = 'Bakover'
end

if strcmp(options.metode,'Bakover')
    Sekant = (FunctionValues(3) - FunctionValues(2))/Timestep;

elseif strcmp(options.metode,'Forover')
    Sekant = (FunctionValues(3) - FunctionValues(2))/Timestep;

elseif strcmp(options.metode,'Senter')
    Sekant = (FunctionValues(3) - FunctionValues(1))/(2*Timestep);

else
    errordlg('Feil metode spesifisert')
    return
end

end
