function loadStruct = getFsmRllCode(rllD)
%GETFSMRLLCODE Load definition of FSM RLL codes from HDD.

%% Create file name
fileName = ['./RllFsmData/fsmD', num2str(rllD), '.mat'];


%% Check if file exists and load
if isfile(fileName)
    loadStruct = load(fileName);
else
    error("Error: The required definition of the FSM RLL code is missing! File: %s ", fileName);
end


end

