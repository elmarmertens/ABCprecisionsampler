function [Xdraw, CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx] = ...
    VARTVPSVprecisionsamplerNaN0const(aaa,invbbb,ccc,y,yNaN,x0,xbar,rndStream,CC,QQ,RR1,arows, acols, asortndx, brows, bcols, bsortndx)
% ABCPRECISIONSAMPLER for case of VAR(p) with missing values and fixed initial conditions
%
% aaa is Ny x Ny x p (x T) (T dimension is optional)
% invbbb is Ny x Ny (x T) (T dimension is optional); invbbb * invbbb' is inverse of variance of VAR residuals
% ccc is measurement vecotr (typically be identity matrix), in general Ny x Ny ( x T)
% x0 is Ny * p vector of initial conditions (p lags of y)
% arguments after rndStream can be empty and will be returned as outputs for use in future calls
% Xdraws is Ny * T vector output (can be shaped to Ny x T)

%% VERSION INFO
% AUTHOR    : Elmar Mertens

% get dimensions
[Ny, T] = size(y);
p       = size(aaa,3);
Nx      = size(x0,1);
Nw      = size(invbbb,2);

if nargin < 9
    CC  = [];
    QQ  = [];
    RR1 = [];
    [arows, acols, asortndx, brows, bcols, bsortndx] = deal([]);
end

if ndims(aaa) == 3
    aaa = repmat(aaa, [1 1 1 T]);
end
if ismatrix(invbbb)
    invbbb = repmat(invbbb, [1 1 T]);
end
if ismatrix(ccc)
    ccc = repmat(ccc, [1 1 T]);
end

NyT   = Ny * T;
NwT   = Nw * T;
NxT   = Nx * T;

%% construct vectorized state space
Y     = reshape(y, NyT, 1);
Ynan  = reshape(yNaN, NyT, 1);
Y     = Y(~Ynan);

%% vectorize input matrices
NxNx         = Nx * Nx;
NxNxT        = NxNx * T;
invbbb       = reshape(invbbb, NxNxT, 1);
ccc          = reshape(ccc, Ny * NxT, 1);

XX0          = repmat(xbar, T, 1);

%% adjust XX0 for initial conditions
% adjust for initial conditions
for k = 1 : p
    xndx                       = (k-1) * Ny + (1 : Ny);
    theseInitialLags           = p - k + 1;
    thisA                      = reshape(aaa(:,:,k:p), Ny, Ny * theseInitialLags);
    thisX0                     = reshape(x0(:,1:p-k+1), Ny * theseInitialLags, 1);
    XX0(xndx)                  = XX0(xndx) + thisA * thisX0; 
end

%% CC and prepare Arows and Brows

if isempty(CC)

    % no pre-allocation of memory here, since to be evaluated only once

    % AA, build sequentially: first unit diagonal
    arows     = 1 : NxT;
    acols     = 1 : NxT;
    % add p lags (sequentially)
    for k = 1 : p
        theserows   = repmat((1 : Nx)', 1 , Nx, T - k);
        theserows   = theserows + permute(Nx * (k : T-1), [1 3 2]);
        arows       = [arows(:); theserows(:)];

        thesecols   = repmat(1 : Nx * (T - k), Nx, 1);
        acols       = [acols(:); thesecols(:)];
    end

    [acols, asortndx] = sort(acols);
    arows             = arows(asortndx);

    % BB
    brows  = repmat((1 : Nx)', 1 , Nw, T);
    brows  = brows + permute(Nw * (0 : T-1), [1 3 2]);
    brows  = reshape(brows, Nx * NwT, 1);
    
    bcols  = repmat(1 : NwT, Nx, 1);
    bcols  = reshape(bcols, NwT * Nx, 1);

    [bcols, bsortndx] = sort(bcols);
    brows             = brows(bsortndx);

    % C
    crows     = repmat((1 : Ny)', 1 , Nx, T);
    crows     = crows + permute(Ny * (0 : T-1), [1 3 2]);
    crows     = crows(:);
    ccols     = repmat(1 : NxT, Ny, 1);
    ccols     = ccols(:);
    CC        = sparse(crows, ccols, ccc, NyT, NxT);
    % drop rows associated with NaN
    CC        = CC(~yNaN,:);
    % perform QR
    [QQ,RR]   = qr(CC');
    [N1, N2]  = size(CC);
    N2        = N2 - N1;
    RR1       = RR(1:N1,1:N1)';

else

    N1        = size(RR1,1);
    N2        = size(QQ,1) - N1;
    
end

QQ1       = QQ(:,1:N1)';
QQ2       = QQ(:,N1+1:end)';

%% sparse builds for BB and AA
% AA
values           = NaN(NxT + NxNx * sum(T - 1 : p),1);
% add unit diagonal element
values(1 : NxT)  = 1;
% add p lags (sequentially)
offset = NxT;
for k = 1 : p
    values(offset + (1 : NxNx * (T-k))) = -reshape(aaa(:,:,k,1+k:T), Nx * Nx * (T - k), 1);
    offset = offset + NxNx * (T-k);
end
values              = values(asortndx);
AA                  = sparse(arows, acols, values, NxT, NxT);

% BB
invbbb  = invbbb(bsortndx);
invBB   = sparse(brows, bcols, invbbb, NxT, NxT);


%% means and innovations
EX        = AA \ XX0;
EY        = CC * EX;

X1tilde   = RR1 \ (Y - EY);

%% precision-based sampler
AAtilde     = invBB * AA;
AAtildeQQ   = AAtilde * QQ;
AAtildeQQ2  = AAtildeQQ(:,N1+1:end);
QQ2invSIG   = AAtildeQQ2' * AAtildeQQ;
invQSIG21   = QQ2invSIG(:,1:N1);
invQSIG22   = QQ2invSIG(:,1+N1:end);


cholinvQSIG22 = chol(invQSIG22, 'lower');


X2hat         = - cholinvQSIG22 \ (invQSIG21 * X1tilde);

Z2draw        = randn(rndStream, N2, 1) + X2hat;
X2draw        = cholinvQSIG22' \ Z2draw;
Xdraw         = EX + QQ1' * X1tilde + QQ2' * X2draw;

