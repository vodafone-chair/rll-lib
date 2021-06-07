function [llrOut] = rllBcjrFct(llrIn, rllD, startState, stopState)
%RLLBCJRFCT This function performs soft-input soft-output demapping of an 
% RLL sequence by implementing the BCJR algorithm.
% ### Inputs ###
% `llrIn` - vector with real-valued input LLRs; one LLR per RLL symbol
% `rllD` - RLL d-constraint; options: rllD={1,2,3,4}
% `startState` - RLL encoder start state; choose `startState=-1` if unknown
% `stopState` - RLL encoder stop state; choose `stopState=-1` if unknown
% ### Outputs ###
% `llrOut` - converted LLRs per bit at the encoder input (soft estimate)
% ### Published in ### 
% P. Neuhaus, M. Dörpinghaus, and G. Fettweis, “Zero-Crossing Modulation for Wideband Systems Employing 1-Bit Quantization and Temporal Oversampling: Implementation and Performance Evaluation,” submitted to IEEE Open Journal of the Communications Society.
% P. Neuhaus, M. Dörpinghaus, H. Halbauer, S. Wesemann, M. Schlüter, F. Gast, and G. Fettweis, “Sub-THz  wideband  system  employing  1-bit  quantization  and  temporal  oversampling,” in Proc. IEEE Int. Conf. Commun. (ICC), Dublin, Ireland, Jun. 2020, pp. 1–7.
% ### References ###
% L. Bahl, J. Cocke, F. Jelinek, and J. Raviv, “Optimal decoding of linearcodes for minimizing symbol error rate (corresp.),”IEEE Trans. Inf.Theory, vol. 20, no. 2, pp. 284–287, Mar. 1974.
% M. Tüchler and A. C. Singer, “Turbo equalization: An overview,” IEEE Trans. Inf. Theory, vol. 57, no. 2, pp. 920–952, Feb. 2011.

%% Load RLL FSM data from HDD
loadStruct = getFsmRllCode(rllD);
outputDef = loadStruct.outputDef;
stateTrans = loadStruct.stateTrans;
Sigma = loadStruct.Sigma;


%% Init
p = log2(size(stateTrans,2));
q = size(outputDef,1);
numStates = size(stateTrans,1);
numOutput = size(llrIn,1)/q;


%% Convert LLR to bit probabilities; Note: LLR = log(P(X=0)/P(X=1))
pIn = ones(size(llrIn,1),2);
pIn(:,1) = 1./(1+exp(-llrIn));
pIn(:,2) = pIn(:,2) - pIn(:,1);

% Reshape for more convenient processing: 
% 1-dim: Probability for -1 and +1
% 2-dim: Vector of n elements corresponding to the length of 1 transition
% 3-dim: Total number of transitions
pIn = reshape(pIn',2,q,numOutput);


%% Obtain transition probabilities Gamma
GammaFull = zeros(numStates,numStates,2^p,numOutput);

% Modified version for confenience
outputDefMod = outputDef;
outputDefMod(outputDef==-1) = 1;
outputDefMod(outputDef==1) = 2;

% Iterate over the number of outputs to compute Gamma
for iter = 1:numOutput
    for oldState = 1:numStates                                              % Iterate over all states
        for curInput = 1:2^p                                                % Iterate over all inputs
            
            % Obtain new state for this transition
            newState = stateTrans(oldState,curInput);

            % Initialize a-priori probability
            tmp = 1/(2^p);

            % Compute probability for current transition by multiplying the
            % probabilities of each bit which would result in this
            % transition
            for qIter = 1:q
                tmp = tmp * pIn(outputDefMod(qIter,oldState,curInput),qIter,iter);
            end
            
            % Update GammaFull
            GammaFull(newState, oldState, curInput, iter) = tmp;

        end
    end
end

% Compute compact version of Gamma
Gamma = squeeze(sum(GammaFull,3));

% Note:
% `GammaFull' corresponds to the state transition probability per input,
% whereas `Gamma' corresponds to the combined state transition probability
% for all inputs.


%% Do BCJR

% Init
alpha = zeros(numStates,numOutput+1);
beta = zeros(numStates,numOutput+1);
app = zeros(numOutput*p,2);

% Define starting and terminating state
if startState == -1
    alpha(:,1) = ones(numStates,1)./numStates;
else
    alpha(startState,1) = 1;
end
if stopState == -1
    beta(:,end) = ones(numStates,1)./numStates;
else
    beta(stopState,end) = 1;
end

for i = 1:numOutput
    % Forward recursion update
    alpha(:,i+1) = Gamma(:,:,i) * alpha(:,i);
    % Normalize
    alpha(:,i+1) = alpha(:,i+1)./sum(alpha(:,i+1));
end

for i = numOutput:-1:1
    % Backward recursion update
    beta(:,i) = Gamma(:,:,i)' * beta(:,i+1);
    % Normalize
    beta(:,i) = beta(:,i)./sum(beta(:,i));
end

for i = 1:numOutput
    for curInput = 1:2^p
        curBin = decToBinFct(curInput-1,p)+1;
        for pIter = 1:p
            curBit = curBin(pIter);
            % APP update
            app((i-1)*p+pIter,curBit) = app((i-1)*p+pIter,curBit) + beta(:,i+1)'*(Sigma(:,:,pIter,curBit).*GammaFull(:,:,curInput,i))*alpha(:,i);
        end
    end
end

% Obtain LLR values
llrOut = log(app(:,1)) - log(app(:,2));


end