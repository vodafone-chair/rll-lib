function [encodedBits, stopState] = rllEncoderFct(inputBits, rllD, startState)
%RLLENCODERFCT Finite state machine based RLL encoder.
% RLL sequence by implementing the BCJR algorithm.
% ### Inputs ###
% `inputBits` - vector bits to be encoded
% `rllD` - RLL d-constraint; options: rllD={1,2,3,4}
% `startState` - RLL encoder start state; default is `startState=1`
% ### Outputs ###
% `encodedBits` - encoded RLL sequence
% `stopState` - final RLL encoder state
% ### Published in ### 
% P. Neuhaus, M. Dörpinghaus, and G. Fettweis, “Zero-Crossing Modulation for Wideband Systems Employing 1-Bit Quantization and Temporal Oversampling: Implementation and Performance Evaluation,” submitted to IEEE Open Journal of the Communications Society.
% P. Neuhaus, M. Dörpinghaus, H. Halbauer, S. Wesemann, M. Schlüter, F. Gast, and G. Fettweis, “Sub-THz  wideband  system  employing  1-bit  quantization  and  temporal  oversampling,” in Proc. IEEE Int. Conf. Commun. (ICC), Dublin, Ireland, Jun. 2020, pp. 1–7.

%% Validate input
assert(all(isreal(inputBits)),'Error: Found imaginary input at the RLL encoder while expexing real binary input!');
assert(all(ismember(inputBits, [0,1])),'Error: Non-binary input at the RLL encoder while expexing real binary input!');

% We are expecting column vectors as input
if size(inputBits,1) == 1 && size(inputBits,2) > 1
    inputBits = transpose(inputBits);
end

% Default option is `startState=1`, if not specified
if nargin < 3 || isempty(startState)
    startState = 1;
end


%% Load RLL FSM data from HDD
loadStruct = getFsmRllCode(rllD);
outputDef = loadStruct.outputDef;
stateTrans = loadStruct.stateTrans;


%% Init
p = log2(size(stateTrans,2));
q = size(outputDef,1);

% Validate input dimension
assert(mod(length(inputBits)/p,1)==0,'Error: RLL encoder input bit vector is of length %d, which cannot be divided into blocks of size p=%d.',length(inputBits),p);

% Allocate memory for output
encodedBits = zeros((q*length(inputBits))/p,1);


%% Encoding
% Initialize encoder to start state
currentState = startState;

% Iteratively encode blocks of `p' bits onto blocks of `q' RLL symbols
for i = 1:length(inputBits)/p
    % Obtain current input bits
    curBits = inputBits((i-1)*p+1:i*p)';
    % Convert to a decimal number which specifies the input
    curBitsDec = binToDecFct(curBits)+1;
    % Obtain the next state based on currentState and current input
    nextState = stateTrans(currentState,curBitsDec);
    % Obtain the ouput based on currentState and current input
    encodedBits((i-1)*q+1:i*q) = outputDef(:,currentState,curBitsDec);
    % Update the next state
    currentState = nextState;   
end

% Save the last state (can be used for decoding)
stopState = currentState;


end
