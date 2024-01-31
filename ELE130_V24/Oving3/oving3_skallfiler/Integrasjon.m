function IntValueNew = Integrasjon(IntValueOld, Timestep, FunctionValues, options)
arguments
    IntValueOld (1,1) double
    Timestep (1,1) double
    FunctionValues (1,2) double
    options.metode (1,:) char = 'Trapes'
end

if strcmp(options.metode,'EulersForover')
    % fyll inn

elseif strcmp(options.metode,'EulersBakover')
    % fyll inn

elseif strcmp(options.metode,'Trapes')
    % fyll inn

else
    errordlg('Feil metode spesifisert')
    return
end

end
