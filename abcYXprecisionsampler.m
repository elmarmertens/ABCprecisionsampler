function [Xdraw, QQ, RR1, Xhat] = abcYXprecisionsampler(Y,XX0,AA,invBB,CC,rndStream,QQ,RR1)
% YXprecisionsampler bare-bones precision sampler for vectorized inputs.
%

%% VERSION INFO
% AUTHOR    : Elmar Mertens

% get dimensions
if nargin < 7
    QQ  = [];
    RR1 = [];
end


%% CC and prepare Arows and Brows

if isempty(QQ)
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

if nargout > 3
    Xhat       = EX + QQ1' * X1tilde + QQ2' * (cholinvQSIG22' \ X2hat);
end
