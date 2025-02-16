function out = approximateBounds(S,phi,opt)
% Approximate bounds of the identified set for different parameters
% via simulation.
% Inputs:
% - restr: structure containing information about restrictions
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options

H = opt.H;
ivar = opt.ivar;
jshock = opt.jshock;
Qdraws = opt.Qdraws;
vma = phi.vma;

n = size(phi.Sigmatr,1);

%% Obtain draws from space of orthonormal matrices satisfying restrictions.
Q0 = drawQs(S,phi,opt);
% Keep relevant column of Q.
Q0 = squeeze(Q0(:,jshock,:));

%% Compute identified sets for impulse responses.
rDraw = zeros(H+1,length(ivar),Qdraws);

for hh = 1:H+1 % For each horizon

    % Extract required rows of horizon-h VMA coefficient matrix.
    Cphi = vma(ivar,:,hh);
    % Multiply by column of Q (for each draw) to obtain impulse responses.
    rDraw(hh,:,:) = Cphi*Q0;

end

% Compute minimum and maximum impulse response over draws of Q.
rmin = min(rDraw,[],3);
rmax = max(rDraw,[],3);

%% Compute identified set for impulse responses to fixed-size shock.
% For each draw of Q, divide impulse responses by impact response of first
% variable.
hDraw = rDraw./rDraw(1,1,:);

% Compute minimum and maximum over draws of Q.
hmin = min(hDraw,[],3);
hmax = max(hDraw,[],3);

%% Compute identified sets for FEVD.
fevhor = zeros(n,1);
fev = zeros(n,H);
conthor = zeros(n,n,Qdraws);
cont = zeros(n,H,Qdraws);

for hh = 1:H % For each horizon
    
    % Extract horizon-h VMA coefficient matrix.
    Cphi = vma(:,:,hh); 
    
    % Compute forecast error variance
    fevhor = fevhor + diag(Cphi*Cphi'); 
    fev(:,hh) = fevhor;
    
    % Compute contribution of first shock to forecast error variance (at
    % each draw of Q).
    cQ = Cphi*Q0;
    for kk = 1:Qdraws
    
        conthor(:,:,kk) = conthor(:,:,kk) + cQ(:,kk)*cQ(:,kk)';
        cont(:,hh,kk) = diag(conthor(:,:,kk));
    
    end
        
end

% Compute draws of FEVD.
fDraw = cont./repmat(fev,[1,1,Qdraws]);

% Compute minimum and maximum of draws over Q.
fmin = permute(min(fDraw,[],3),[2,1]);
fmax = permute(max(fDraw,[],3),[2,1]);

%% Compute identified set for covariance between proxy and first shock.
dphi = phi.d'*phi.Sigmatr;
gDraw = dphi*Q0;
gmin = min(gDraw);
gmax = max(gDraw);

%% Collect outputs.
out.rmin = rmin;
out.rmax = rmax;
out.hmin = hmin;
out.hmax = hmax;
out.fmin = fmin;
out.fmax = fmax;
out.gmin = gmin;
out.gmax = gmax;

end

