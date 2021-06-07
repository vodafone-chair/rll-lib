function R = getFsmRllAutocorrFct(rllD, lenAutocorr, type)
%GETFSMRLLAUTOCORRFCT This function computes the autocorrelation function
% of the FSM RLL codes.
% ### Inputs ###
% `rllD` - RLL d-constraint; options: rllD={1,2,3,4}
% `lenAutocorr` - length of the computed autocorrelation function
% `type` - type of the computed autocorrelation; options: type={'one-sided', 'two-sided'}
% ### Outputs ###
% `R` - autocorrelation of the FSM RLL codes
% ### Published in ###
% P. Neuhaus, M. Dörpinghaus, and G. Fettweis, “Zero-Crossing Modulation for Wideband Systems Employing 1-Bit Quantization and Temporal Oversampling: Implementation and Performance Evaluation,” submitted to IEEE Open Journal of the Communications Society.
% P. Neuhaus, M. Dörpinghaus, H. Halbauer, V. Braun, and G. Fettweis, “On the spectral efficiency of oversampled 1-bit quantized systems for wideband LOS channels,” in Proc. IEEE Int. Symp. on Personal, Indoor and Mobile Radio Commun. (PIMRC), London, U.K., Aug. 2020. 


%% Validate inputs
if nargin < 3 || isempty(type)
    type = 'two-sided';
end


%% Load RLL FSM data from HDD
loadStruct = getFsmRllCode(rllD);
Gamma = loadStruct.Gamma;
P = loadStruct.P;


%% Init
q = size(Gamma,2);


%% Compute the stationary distribution
% Find right eigenvector corresponding to eigenvalue lambda=1
[~,D,W] = eig(P);
% Allow for numerical inaccurary
idx = abs(diag(D)-1)<1e-9;
e = W(:,idx);
% Normalize total probability to 1 (unit norm eigenvector)
e = e./sum(e);

% Zero-padding is implemented implicitly
curP = eye(length(P));

% Preallocate
R_k = zeros(size(Gamma,2),size(Gamma,2),lenAutocorr);

% Compute block-wise autocorrelation matrix
for m = 1:ceil(lenAutocorr/q)
    R_k(:,:,m) = Gamma'*diag(e)*curP*Gamma;
    curP = P * curP;
end


%% Phase averaging
for k = 1:ceil(lenAutocorr/q)
    for l = 1:q
        R((k-1)*q+l) = 0;
        for i = 1:q-(l-1)
            R((k-1)*q+l) = R((k-1)*q+l) + R_k(i,(l-1)+i,k);
        end
        for i = q-(l-1)+1:q
            R((k-1)*q+l) = R((k-1)*q+l) + R_k(i,(l-1)+i-q,k+1);
        end
        % Divide by q
        R((k-1)*q+l) = R((k-1)*q+l) / q;
    end
end


%% Convert to desired output format
switch type
    case 'one-sided'
        R = R(1:lenAutocorr);
    case 'two-sided'
        if mod(lenAutocorr,2) == 0
            lenAutocorr = lenAutocorr + 1;
            warning('Two-sided autocorrelations are always computed with odd length. Hence, the output length has been increased to %d', lenAutocorr);
        end
        R = [fliplr(R(2:1+(lenAutocorr-1)/2)),R(1:1+(lenAutocorr-1)/2)];
    otherwise
        error('Error: Undefined type=`%s`. Expected `one-sided` or `two-sided`.',type);
end


end

