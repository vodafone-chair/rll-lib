%% Clean up
close all;
clear all;
clc;


%% Init
snr_dB = 5;                                                                 % define SNR in dB
rllD = 4;                                                                   % RLL constraint `d' (options: 1/2/3/4)
seqLen = 1e4+2;                                                             % binary input sequence length
rng('shuffle');                                                             % init random number generator
inputBits = randi([0 1],seqLen,1);                                          % generate random input bits


%% Encode bits to RLL sequence
startState = 1;                                                             % State in which the RLL encoder starts, typically `1' is a good choice
[encodedBits, stopState] = ...                                              % RLL encoding
    rllEncoderFct(inputBits, rllD, startState);                             


%% AWGN channel
% Compute noise variance
noiseVar = 10^(-snr_dB/10);
% Generate noise
n = noiseVar .* randn(length(encodedBits),1);
% Add noise
y = encodedBits + n;
% Convert observations to LLRs
llrIn = - 2/noiseVar .* y;


%% RLL decoding
llrOut = rllBcjrFct(llrIn, rllD, startState, stopState);                    % Soft RLL decoding (BCJR)
bitsHatHard = rllViterbiFct(y, rllD, startState, stopState);                % Hard RLL decoding (Viterbi)

% Hard demapping of the soft-output from BCJR
bitsHatSoft = ones(seqLen,1);
bitsHatSoft(llrOut>0) = 0;


%% Eval uncoded BER
uncodedBER(1) = numel(find(bitsHatHard~=inputBits))/seqLen;
uncodedBER(2) = numel(find(bitsHatSoft~=inputBits))/seqLen