function mickey_mouse_path = importMickey(filename, dataLines)
%IMPORTFILE2 Import data from a text file
%  MICKEY_MOUSE_PATH = IMPORTFILE2(FILENAME) reads data from text file
%  FILENAME for the default selection. Returns the numeric data as a Nx2 array.
%
%  MICKEY_MOUSE_PATH = IMPORTFILE2(FILENAME, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME.

%% Input handling
if nargin < 2
    dataLines = [1, 230];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

opts.DataLines = dataLines;
opts.Delimiter = ";";

opts.VariableNames = ["x", "y"];
opts.VariableTypes = ["double", "double"];

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

opts = setvaropts(opts, ["x", "y"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["x", "y"], "DecimalSeparator", ",");
opts = setvaropts(opts, ["x", "y"], "ThousandsSeparator", ".");

%% Import the data
T = readtable(filename, opts);

%% Convert to output type (matrice 230x2)
mickey_mouse_path = table2array(T);

end