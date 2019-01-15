function [ P ] = importpowermeter( filename )
    %IMPORTPOWERMETER Import data from a Thorlabs PM100D power meter.
    %   Imports data in the format saved to the memory card of a Thorlabs
    %   PM100D power meter. Returns data and metadata in a structure.
    %
    % Iain Robinson, School of Geosciences, University of Edinburgh
    % 2015-06-04
    
    % Check that the file exists.
    if ~exist(filename, 'file')
        error('The file "%s" was not found.', filename);
    end
    
    % Read the file.
    text = fileread(filename);

    % Parse the metadata using regular expressions.
    P.sensorSerial = char(regexp(text, '^(\S*)\s', 'tokens', 'once'));
    datetimeString = char(regexp(text, '(\d\d\d\d-\d\d-\d\d \d\d:\d\d:\d\d)', 'tokens', 'once'));
    % Convert the datetime string format from the file's format to
    % MATLAB's.
    P.datetime = datestr(datenum(datetimeString));
    P.valueUnit = char(regexp(text, 'value unit \[(\S+)\]', 'tokens', 'once'));
    P.timeUnit = char(regexp(text, 'time unit \[(\S+)\]', 'tokens', 'once'));
    P.wavelength = str2double(regexp(text, 'wavelength\s+([\d\.]+)', 'tokens', 'once'));
    P.wavelengthUnit = char(regexp(text, 'wavelength\s+[\d\.]+(\S+)', 'tokens', 'once'));
    P.range = char(regexp(text, 'range\s+(\S+)', 'tokens', 'once'));

    % The data region of the text file extends from the first line beginning
    % with a number, plus sign or minus sign, to the end of the file.
    dataRegion = char(regexp(text, '^[-\+\d].*', 'match', 'lineanchors'));
    if isempty(dataRegion)
        error('Could not find data section of file %s.', filename);
    end

    % Scan with textscan.
    scannedData = textscan(dataRegion, '%f %f');

    % Copy into variables.
    P.time = scannedData{2};
    P.data = scannedData{1};
end

