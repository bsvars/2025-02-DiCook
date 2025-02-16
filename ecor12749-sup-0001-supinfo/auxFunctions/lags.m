function YLag = lags(Y,p)
% Function generates p lags of matrix Y.
% Inputs:
% Y: matrix
% p: vector containing lag orders (p >= 0).

nLags = length(p); % Number of lags to generate

[nObs,nVars] = size(Y); 

YLag = NaN(nObs,nLags*nVars); % Storage

for ii = 1:nLags
    
    Lag = p(ii);
    
    cols = (nVars*(ii-1)+1):(ii*nVars); % Column indices
    
    YLag((Lag+1):end,cols) = Y(1:(end-Lag),:);
      
end

end
