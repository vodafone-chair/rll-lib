function bin = decToBinFct(dec, n)
%DECTOBINFCT Convert decimal input to binary output of length `n` with right
% MSB.
% ### Reference ###
% Parts are taken from https://de.mathworks.com/matlabcentral/fileexchange/26447-efficient-convertors-between-binary-and-decimal-numbers

%% Init
bin = zeros(n,1);


%% Zero Handler
if dec==0
    return
end


%% Conversion
c = floor(log2(dec)) + 1;
for i = c:-1:1
    r = floor(dec / 2);
    bin(n-(c-i)) = dec - 2*r;
    dec = r;
end


end

