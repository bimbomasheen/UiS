function[FilteredValue] = IIR_filter(OldFilteredValue,Measurement,Parameter)

FilteredValue = (1-Parameter) * OldFilteredValue + Parameter * Measurement;

end