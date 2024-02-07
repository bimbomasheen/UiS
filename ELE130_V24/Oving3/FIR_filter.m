function[FilteredValue] = FIR_filter(Measurement,NoMeas)

persistent Measurements

% Initialize Measurements if it's empty
if isempty(Measurements)
    Measurements = Measurement;
else
    Measurements = [Measurements, Measurement]; % Append the new measurement
end

% Adjust NoOfMeasurements if there are not enough measurements
if length(Measurements) < NoMeas
    NoMeas = length(Measurements);
end

% Compute the filtered value as the average of the last NoOfMeasurements measurements
start_index = max(1, length(Measurements)-NoMeas+1);
FilteredValue = sum(Measurements(start_index:end)) / NoMeas;

end

