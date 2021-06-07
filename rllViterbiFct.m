function [bitsHat] = rllViterbiFct(y, rllD, startState, stopState)
%RLLVITERBIFCT This function performs Viterbi demapping of an 
% RLL sequence using the Hamming distance criterion.
% ### Inputs ###
% `y` - vector with real observations; one observation per RLL symbol
% `rllD` - RLL d-constraint; options: rllD={1,2,3,4}
% `startState` - RLL encoder start state; choose `startState=-1` if unknown
% `stopState` - RLL encoder stop state; choose `stopState=-1` if unknown
% ### Outputs ###
% `bitsHat` - estimated bits at the encoder input (hard estimate)
% ### Published in ###
% P. Neuhaus, D. M. V. Melo, L. T. N. Landau, R. C. de Lamare, and G. Fettweis. “Zero-Crossing Modulations for a Multi-User MIMO Downlink with 1-Bit Temporal Oversampling ADCs”. In Proc. European Signal Proc. Conf. (EUSIPCO). Dublin, Ireland, Sept. 2021.
% ### References ###
% A. Viterbi, “Error bounds for convolutional codes and an asymptotically optimum decoding algorithm,” IEEE Trans. Inf. Theory, vol. 13, no. 2, pp. 260–269, Apr. 1967.
% G. D. Forney, “The Viterbi algorithm,” Proc. IEEE, vol. 61, no. 3, pp. 268–278, Mar. 1973.


%% Validate inputs
if ~all(abs(y)==1)
    % If input does not consist of 1-bit measurements, take sign only due
    % to minimum Hamming distance decoding
    y = sign(y);
end

% Check if input is only real valued
assert(all(isreal(y)),'Error: Found imaginary input at the RLL decoder input while expexing only real inputs!');

% We are expecting column vectors as input, otherwise transpose
if size(y,1) == 1 && size(y,2) > 1
    y = transpose(y);
end


%% Load RLL FSM data from HDD
loadStruct = getFsmRllCode(rllD);
outputDef = loadStruct.outputDef;
stateTrans = loadStruct.stateTrans;


%% Init
p = log2(size(stateTrans,2));
q = size(outputDef,1);
numTransitions = size(y,1)/q;
numStates = size(stateTrans,1);

survivorSeq = Inf .* ones(numTransitions,numStates);

% Init start state
if startState ~= -1
    pathMetric = Inf .* ones(numStates,1);
    pathMetric(startState) = 0;
else
    pathMetric = zeros(numStates,1);
end


%% Viterbi decoding
for iter = 1:numTransitions
    
    % Get current observation
    curObs = y(1+(iter-1)*q:iter*q);
    
    % Compute current branch metrics
    branchMetric = Inf .* ones(numStates,numStates,2^p);
    
    % Iterate over all states
    for curState = 1:numStates
        
        % Iterate over all decimal inputs
        for curInputDec = 1:2^p
            
            % Obtain the next state based on current state and current input
            nextState = stateTrans(curState,curInputDec);

            % Compute branch metric (Hamming distance) for this transition
            branchMetric(nextState,curState,curInputDec) = pathMetric(curState) + sum(abs(curObs - outputDef(:,curState,curInputDec)))/2;
            
        end
        
    end
    
    % Min search over inputs
    [minBranchMetric, decIdx] = min(branchMetric,[],3);
    
    % Update path metric and find survivor state indices by min search
    % over current states
    [pathMetric, survivorIdx] = min(minBranchMetric,[],2);
    
    % Obtain corresponding decimal inputs
    myIdx = sub2ind(size(decIdx),1:size(decIdx,1),survivorIdx');
    decVec = decIdx(myIdx)';    
    
    % Store survivors sequences
    survivorSeqOld = survivorSeq;
    for i = 1:numStates
        if isfinite(pathMetric(i))
            decOut = decVec(i);
            if iter > 1
                survivorSeq(1:iter*p,i) = [survivorSeqOld(1:(iter-1)*p,survivorIdx(i));decToBinFct(decOut-1,p)];
            else
                survivorSeq(1:iter*p,i) = decToBinFct(decOut-1,p);
            end
        end
    end
    
end


%% Choose survivor sequence
if stopState ~= -1                                                          % If a terminal state has been specified, choose sequence correspondingly
    bitsHat = survivorSeq(:,stopState);
else                                                                        % Otherwise select stop state as state with minimum path metric
    [~,idxMin] = min(pathMetric);
    bitsHat = survivorSeq(:,idxMin);
end


end