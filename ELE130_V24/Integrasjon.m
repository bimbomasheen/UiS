function IntValueNew = Integrasjon(IntValueOld, Timestep, FunctionValues, options)
arguments
    IntValueOld (1,1) double
    Timestep (1,1) double
    FunctionValues (1,2) double
    options.metode (1,:) char = 'Trapes'
end

if strcmp(options.metode,'EulersForover')
    IntValueNew = IntValueOld + Timestep * FunctionValues(2);

elseif strcmp(options.metode,'EulersBakover')
    IntValueNew = IntValueOld + Timestep * FunctionValues(1);

elseif strcmp(options.metode,'Trapes')
    IntValueNew = IntValueOld + Timestep * 0.5 * (FunctionValues(1) + FunctionValues(2));

else
    errordlg('Feil metode spesifisert')
    return
end

end
