function dec = binToDecFct( bin )
%BINTODECFCT Converts a binary input vector to a decimal number.
% ### Reference ###
% Parts are taken from https://de.mathworks.com/matlabcentral/fileexchange/26447-efficient-convertors-between-binary-and-decimal-numbers

%% Conversion
dec = sum(bin.*(2.^(size(bin,2)-1:-1:0)),2);