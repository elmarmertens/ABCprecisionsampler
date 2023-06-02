function [YbarDraw, YgapDraw, Abar, agaprows, agapcols, agapsortndx, brows, bcols] = ...
    trendcyclePrecisisionsampler(y, agap, invbgap, invbbar, ybar0, rndStream, Abar, agaprows, agapcols, agapsortndx, brows, bcols)

%% VERSION INFO
% AUTHOR    : Elmar Mertens

if nargin < 7
    [Abar, agaprows, agapcols, agapsortndx, brows, bcols] = deal([]);
end
% get dimensions
[Ny, T] = size(y);
p       = size(agap,3);

if ndims(agap) <= 3
    agap = repmat(agap, [1 1 1 T]);
end
if ismatrix(invbgap)
    invbgap = repmat(invbgap, [1 1 T]);
end
if ismatrix(invbbar)
    invbbar = repmat(invbbar, [1 1 T]);
end

NyT = Ny * T;

%% construct vectorized state space
Y      = reshape(y, NyT, 1);
Ybar0  = sparse(1:Ny, 1, ybar0, NyT, 1);

NyNy   = Ny * Ny;
Nagap  = NyT + NyNy * p * (T - p) + sum(NyNy * (1 : p - 1));
Nb    = NyNy * T;
    
if isempty(Abar)
    % Abar
    Abar = spdiags([-1 * ones(NyT,1), ones(NyT,1)], [-Ny, 0], NyT,NyT);

    % Agap
    [agaprows, agapcols] = deal(NaN(Nagap, 1));
    agaprows(1:NyT) = 1:NyT;
    agapcols(1:NyT) = 1:NyT;
    offset = NyT;
    for k = 1 : p
        theserows   = repmat((1 : Ny)', 1 , Ny, T - k);
        theserows   = theserows + permute(Ny * (k : T-1), [1 3 2]);

        thesecols  = repmat(1 : Ny * (T - k), Ny, 1);

        agaprows(offset + (1 : NyNy * (T-k))) = theserows(:);
        agapcols(offset + (1 : NyNy * (T-k))) = thesecols(:);

        offset = offset + NyNy * (T-k);
    end
    % sort Agap indices
    ndx = sub2ind([NyT, NyT], agaprows, agapcols);
    [~, agapsortndx] = sort(ndx);
    agaprows         = agaprows(agapsortndx);
    agapcols         = agapcols(agapsortndx);

    % Bbar and Bgap use same row indices
    brows = repmat((1:Ny)', 1, Ny) + permute(Ny * (0 : T - 1), [1 3 2]);
    bcols = repmat(1:Ny, Ny, 1) + permute(Ny * (0 : T - 1), [1 3 2]);
    brows = reshape(brows, Nb, 1);
    bcols = reshape(bcols, Nb, 1);
    % no sorting needed
    % ndx   = sub2ind([NyT, NyT], brows, bcols);
    % [~, bsortndx] = sort(ndx);
    % brows = brows(bsortndx);
    % bcols = bcols(bsortndx);

end

avalues           = ones(Nagap,1);
offset = NyT;
for k = 1 : p
    avalues(offset + (1 : NyNy * (T-k))) = -reshape(agap(:,:,k,1+k:T), NyNy * (T - k), 1);
    offset = offset + NyNy * (T-k);
end
avalues             = avalues(agapsortndx);
Agap                = sparse(agaprows, agapcols, avalues, NyT, NyT);

invBbar = sparse(brows, bcols, invbbar(:), NyT, NyT);
invBgap = sparse(brows, bcols, invbgap(:), NyT, NyT);

Abar      = invBbar * Abar;
P0bar     = Abar' * Abar;

Agap      = invBgap * Agap;
P0tilde   = Agap' * Agap;

Ybar0     = invBbar * Ybar0; % Ytilde0 is assumed zero

Pbar         = P0bar + P0tilde;
sqrtPbar     = chol(Pbar, 'lower');

Pbarmu       = Abar' * Ybar0 + P0tilde * Y;
sqrtPmu      = sqrtPbar \ Pbarmu;
Zdraw        = randn(rndStream, NyT, 1);
YbarDraw     = transpose(sqrtPbar) \ (sqrtPmu + Zdraw);
if nargout > 1
    YgapDraw     = Y - YbarDraw;
end