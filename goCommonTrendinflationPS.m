%% trend inflation model with state space sampling via Precision-based sampling


%% clear workspace
clear
close all
fclose all;
clc

rng(061222); %  fix random seed


%#ok<*NASGU>
%#ok<*UNRCH>
%#ok<*DATNM>
%#ok<*DATST>

%% load toolboxes
path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/


%% defaults
samStart  = [];
samStop   = [];

datalabel     = 'INFTRMSRV';
np            = 12;
PAIpriortheta = [0.2^2 0.5^2 2 100];
PAIpriortheta = [0.2^2 0.5^2 4 100];
p             = np;
NfcstHorizons = np * 10;


fontsize= 16;

% MCMC
Nmcmc     = 1e4;
burnin    = Nmcmc;

doSingleThread = false;

%% SED-PARAMETERS-HERE

% this cell can be used to patch in alternative parameter choices via an "sed" script
% note: only parameters with defaults set prior to this cell can be altered


%% process settings
maxNumCompThreads('automatic');
if doSingleThread
    % enforce single threaded computations
    usedThreads = 1;
    availableThreads = maxNumCompThreads(usedThreads);
    fprintf('Using 1 of %d available threads.\n', availableThreads)
else
    usedThreads =  maxNumCompThreads('automatic');
    availableThreads = usedThreads;
    fprintf('Using all %d available threads.\n', usedThreads)
end


tailquants = normcdf([-1 1]) * 100;


%% load data
FREDtable = readtable(sprintf('%s.csv', datalabel));

Yinflation  = table2array(FREDtable(:,2:end));
dates       = datetime(FREDtable.dates, 'ConvertFrom', 'datenum');
Ycode       = FREDtable.Properties.VariableNames(2:end);

T     = length(dates);
N     = size(Yinflation,2);

rndStream = getDefaultStream;

%% initwrap
switch np
    case 4
        datefmt = 'yyyyQQ';
    case 12
        datefmt = 'yyyymm';
    otherwise
        error('np == %d not yet supported', np)
end

titlename = sprintf('commontrendinflation-%s-Ndraws%s', datalabel, erase(sprintf('%1.0e', Nmcmc), '+'));
if ~isempty(samStart)
    titlename = sprintf('%s-%s', titlename, datestr(samStart, datefmt));
end
if ~isempty(samStop)
    titlename = sprintf('%s-%s', titlename, datestr(samStop, datefmt));
end

execLabel = 'PS';

if doSingleThread
    execLabel = strcat(execLabel, 'singlethreaded');
else
    execLabel = strcat(execLabel, 'multithreaded');
end
titlename = strcat(titlename, '-', execLabel);
initwrap

if isempty(wrap) && (~isdesktop || ispc)
    initwrap
end

%% cut sample
if isempty(samStart)
    samStart = dates(1);
end
if isempty(samStop)
    samStop = dates(end);
end
sam    = (dates >= samStart) & (dates <= samStop);

dates         = dates(sam);
Yinflation    = Yinflation(sam,:);
T             = length(dates);


%% plot data
% thisfig = figure;
% hold on
% for n = 1 : N
%     if iscompact(Yinflation(:,n))
%         plot(dates, Yinflation(:,n), '-')
%     else
%         plot(dates, Yinflation(:,n), 'd')
%     end
% end
% xtickdates(dates)
% wrapthisfigure(thisfig, sprintf('alldata-%s', datalabel), wrap);
%
% thisfig = figure;
% hold on
% plot(dates, meanK(Yinflation,12,1))
% plot(dates, mean(meanK(Yinflation,12,1),2, 'omitnan'), 'k-d')
% xtickdates(dates)
% wrapthisfigure(thisfig, sprintf('alldataMA12-%s', datalabel), wrap);

%% prepare state space
Nbar    = 1;
Ngap    = N;
Nstates = Nbar + Ngap * p;
Nshocks = Nbar + Ngap;

ndxbar            = 1 : Nbar;

ndxgap            = Nbar + (1 : Ngap);
ndxgapstates      = Nbar + 1 : Nstates;
ndxgapshocks      = Nbar + (1:Ngap);


%% prepare state space inputs for precision-sampling
aaa0                  = zeros(Nbar + Ngap, Nbar + Ngap, p); % note: 3rd dim of aaa will later be sorted from p to 1
aaa0(ndxbar,ndxbar,1) = eye(Nbar); % random walk

invbbb               = zeros(Nbar + Ngap,Nbar + Ngap,T);
ccc                  = [ones(N,1) eye(N)];

X00                  = sparse((Nbar + Ngap) * p, 1);
invsqrtSigma00       = [1e-4 * ones(Nbar,1); 1e-1 * ones(Ngap,1)];
invsqrtSigma00       = repmat(invsqrtSigma00, p, 1);
invsqrtSigma00       = spdiags(invsqrtSigma00, 0, (Nbar + Ngap) * p, (Nbar + Ngap) * p);


%%  trend SV prior
ElogSVbar0          = 0;
VlogSVbar0          = 100;
varlogSVbar0        = 1;
varlogSVbarDof      = 3;
varlogSVbarT        = (varlogSVbarDof - 2) .* varlogSVbar0;

[KSC, KSCt, logy2offset] = getKSC10values(T, Nbar); % using 10-point grid as in Omori et al (2007, JoE)

%% gap variance prior
OMEGA0gap       = eye(Ngap);
dof0gap         = 2 + Ngap;
OMEGA0gapT      = (dof0gap - 2) * OMEGA0gap;

%% VAR prior
% the prior is for transpose(PAI) when PAI is N times N*p
Nregressors  = Ngap * p + 1;
Npai         = Ngap * Nregressors;
PAIpriormean = zeros(Ngap * Nregressors, 1);
PAIpriorvar  = NaN(Npai,1); % only diagonal elements are stored
ndx = 0;
for i = 1 : Ngap
    for l = 1 : p
        for j = 1 : Ngap
            ndx  = ndx + 1;
            if (i==j)
                PAIpriorvar(ndx) = PAIpriortheta(1) / (l^PAIpriortheta(3));
            else
                PAIpriorvar(ndx) = PAIpriortheta(1) * PAIpriortheta(2)  / (l^PAIpriortheta(3));
            end
        end % j
    end % l
    % const
    ndx  = ndx + 1;
    PAIpriorvar(ndx) = PAIpriortheta(4);
end % i
invPAIpriorvar             = diag(1 ./ PAIpriorvar);
invPAIpriorvarPAIpriormean = invPAIpriorvar * PAIpriormean(:);



%% prepare MCMC

% allocate memory for draws
drawPAI           = NaN(Ngap, Ngap * p, Nmcmc);
drawMaxRoot       = NaN(1, Nmcmc);
drawSQRTOMEGAgap  = NaN(Ngap, Ngap, Nmcmc);
drawStateGap      = NaN(T, Ngap, Nmcmc);
drawStateBar      = NaN(T, Nbar, Nmcmc);
StateGap          = NaN(T+p, Ngap);

drawY             = NaN(T, N); % to track missing data

drawlogSVbar      = NaN(Nbar, T, Nmcmc);
drawlogvarSVbar   = NaN(Nbar, Nmcmc);

drawXjumpoff      = NaN(Nbar+Ngap*p,Nmcmc); % jumpoff vector for forecast sims


% gap VAR mean
drawPAI0          = NaN(Ngap, Nmcmc);
drawY0            = NaN(1,Ngap, Nmcmc); % dimensions chosen to match drawStateBar
% prepare constants
iota         = ones(N,1);

Ingap        = eye(Ngap);
unity        = 1 - 1e-4;
Ntrials      = 1e2;

% prepare VAR companion
kompanion                    = zeros(N * p);
kompanion(N+1:end,1:N*(p-1)) = eye(N*(p-1));

% init Omega and PAI
sqrtOMEGAgap  = Ingap;
logSV         = zeros(Nbar,T);
logSV0        = 0;
varlogSVbar   = varlogSVbar0;

PAI           = zeros(N, N * p);
Y0            = zeros(Ngap,1);
Yinflation    = transpose(Yinflation);

%% prepare handling of missing obs
Ynan             = isnan(Yinflation);
Yinflation(Ynan) = 0;


%% pre-construct precision-sampling objects
tic % startmeasuring execution time

aaa                                = aaa0;
PAIdummy                           = .1 * ones(N, N * p); % dummy value, needs to be different from zero for proper index construction
aaa(Nbar+(1:Ngap),Nbar+(1:Ngap),:) = reshape(PAIdummy, Ngap, Ngap, p);
aaa                                = aaa(:,:,p:-1:1);

invbbb(ndxbar,ndxbar,:)       = eye(Nbar) .* permute(exp(-.5 * logSV), [1 3 2]);
invsqrtOMEGAgap               = randn(rndStream,Ngap,Ngap); % just for init
invbbb(ndxgap,ndxgapshocks,:) = repmat(invsqrtOMEGAgap, [1 1 T]);
[~, CCC, QQQ, RRR1, arows, acols, a0ndx, asortndx, brows, bcols, b0ndx, bsortndx] = ALBCprecisionsamplerNaN(aaa,invbbb,ccc,Yinflation,Ynan,X00,invsqrtSigma00,rndStream);


%% loop over MCMC steps

progressbar(0);
for thisMCMCstep = -burnin : Nmcmc


    %% filter data
    aaa                                = aaa0;
    aaa(Nbar+(1:Ngap),Nbar+(1:Ngap),:) = reshape(PAI, Ngap, Ngap, p);
    aaa                                = aaa(:,:,p:-1:1);

    invbbb(ndxbar,ndxbar,:)       = eye(Nbar) .* permute(exp(-.5 * logSV), [1 3 2]); % w/Nbar=1 could also do "1 ./ exp(.5 * logSV)"
    invsqrtOMEGAgap               = Ingap / sqrtOMEGAgap;
    invbbb(ndxgap,ndxgapshocks,:) = repmat(invsqrtOMEGAgap, [1 1 T]);

    Ytilde = Yinflation - Y0;

    Xdraw = ALBCprecisionsamplerNaN(aaa,invbbb,ccc,Ytilde,Ynan,X00,invsqrtSigma00,rndStream,...
        CCC,QQQ,RRR1,arows,acols,a0ndx,asortndx,brows,bcols,b0ndx,bsortndx);


    States                = reshape(Xdraw, Nbar + Ngap, T + p);
    StateBar              = States(ndxbar,p+(1:T));
    shockBar              = diff(States(ndxbar,p+(0:T)), [], 2);
    StateGap              = transpose(Y0 + States(ndxgap, :));

    gapjumpoff              = fliplr(States(ndxgap,T-p+1:T)); % first collect gaps
    gapjumpoff              = reshape(gapjumpoff, Ngap * p, 1);
    Xjumpoff                = cat(1, States(ndxbar,T+p), gapjumpoff);

    ydraw                   = Y0 + ccc * States(:,p+(1:T));

    %% estimate VAR and draw coeff
    [X, Y]           = lag4VAR(StateGap, p, true);
    iSigmaResid      = invsqrtOMEGAgap' * invsqrtOMEGAgap;

    % generate stable draw
    isstabledraw = false;
    while ~isstabledraw

        [thisPAI, shockGap]  = bayesVectorRegressionGibbsDraw1(Y, X, iSigmaResid, invPAIpriorvarPAIpriormean, invPAIpriorvar, rndStream);
        % note: thisPAI is Nregressors times N

        PAI       = transpose(thisPAI(1:end-1,:));
        PAI0      = transpose(thisPAI(end,:));

        kompanion(1:N,:) = PAI;
        abseigenvalues   = abs(eig(kompanion));
        isstabledraw     = all(abseigenvalues < 1);

    end

    Y0    = (Ingap - sum(reshape(PAI,Ngap,Ngap,p),3)) \ PAI0;

    %% gap shock variances
    sqrtOMEGAgap = bayesSQRTVCVgibbsDraw1(OMEGA0gapT, dof0gap, shockGap, rndStream);

    %% SV shock variance
    logSVshock          = diff([logSV0, logSV], [], 2);
    varlogSVbar         = igamVarianceDraw(logSVshock, varlogSVbarT, varlogSVbarDof);

    %% trend shock SV
    logy2               = log(shockBar.^2 + logy2offset);
    [logSV, logSV0]     = StochVolKSC(logy2, logSV, sqrt(varlogSVbar), ElogSVbar0, VlogSVbar0, KSC, KSCt, Nbar, T, rndStream);

    %% collect draws
    if thisMCMCstep > 0
        drawPAI(:,:,thisMCMCstep)           = PAI;
        drawPAI0(:,thisMCMCstep)            = PAI0;
        drawY0(1,:,thisMCMCstep)            = Y0';
        drawMaxRoot(:,thisMCMCstep)         = max(abseigenvalues);
        drawSQRTOMEGAgap(:,:,thisMCMCstep)  = sqrtOMEGAgap;
        drawStateBar(:,:,thisMCMCstep)      = transpose(StateBar);
        drawStateGap(:,:,thisMCMCstep)      = StateGap(p+1:end,:);
        drawXjumpoff(:,thisMCMCstep)        = Xjumpoff;

        drawY(:,:,thisMCMCstep)             = transpose(ydraw);

        drawlogSVbar(:,:,thisMCMCstep)      = logSV;
        drawlogvarSVbar(:, thisMCMCstep)    = varlogSVbar;
    end

    progressbar((thisMCMCstep + burnin + 1) / (burnin + Nmcmc + 1))

end
MCMCtime = toc;
fprintf('MCMC took %8.2f (with %s)\n', MCMCtime, execLabel)

%% clear MCMC helper variables
clear Y X kompanion
Yinflation(Ynan) = NaN;
Yinflation       = transpose(Yinflation);

checkdiff(Yinflation - drawY); % ignores NaN


%% report maxroot of VAR
thisfig = figure;
histogram(drawMaxRoot, 100)
wrapthisfigure(thisfig, sprintf('VARmaxroot-%s', datalabel), wrap);

%% trend variance
drawSV = exp(.5 * drawlogSVbar);
medSV  = mean(drawSV, 3);
tailSV = prctile(drawSV, normcdf([-1 1]) * 100, 3);

for n = 1 : Nbar
    thisfig = figure;
    hold on
    plot(dates, medSV(n,:), 'k-', 'LineWidth', 2);
    plot(dates, squeeze(tailSV(n,:,:)), 'k-', 'LineWidth', 1);
    ylim0
    wrapthisfigure(thisfig, sprintf('SVbar%d-%s',  n, datalabel), wrap);
end

%% report StateBar
medStateBar  = median(drawStateBar, 3);
tailStateBar = permute(prctile(drawStateBar, normcdf([-1 1]) * 100, 3), [1 3 2]);

for n = 1 : Nbar
    thisfig = figure;
    hold on
    h = plot(dates, medStateBar(:,n), '-', 'color', [0 0 0], 'linewidth', 2);
    plot(dates, tailStateBar(:,:,n), '-.', 'color', [0 0 0], 'linewidth', 1)
    xtickdates(dates)
    ylim0
    grid on
    wrapthisfigure(thisfig, sprintf('Trend%d-%s',  n, datalabel), wrap);
end

%% report intercepts
thesedraws = squeeze(drawY0);
medY0      = median(thesedraws, 2);
thisfig = figure;
hold on
bar(1:N, medY0)
xticks(1:N)
xticklabels(Ycode)
grid on
wrapthisfigure(thisfig, sprintf('Y0-%s',  datalabel), wrap);

%% report Trends
drawTrends    = drawStateBar + drawY0;

medTrends     = median(drawTrends, 3);
tailTrends    = prctile(drawTrends, tailquants, 3);

medY          = median(drawY,3);
tailY         = prctile(drawY, normcdf([-1 1]) * 100, 3);

for n = 1 : N
    thisfig = figure;
    hold on
    plotCI(medY(:,n), squeeze(tailY(:,n,:)), dates, [], ':', 'color', [.25 .25 .25]);
    h = plot(dates,  medTrends(:,n), '-', 'linewidth', 2);
    if iscompact(Yinflation(:,n))
        hData = plot(dates,Yinflation(:,n), ':', 'color', 'black', 'linewidth', 1);
    else
        hData = plot(dates,Yinflation(:,n), 'd', 'color', 'black', 'linewidth', 1);
    end
    hDataMA = plot(dates, meanK(Yinflation(:,n), np, 1), '-', 'color', 'black', 'linewidth', 1);
    legend([h hData hDataMA], 'Trend', 'Data', sprintf('Data MA(%d)', np))
    xtickdates(dates)
    ylim0
    grid on
    title(Ycode{n})
    wrapthisfigure(thisfig, sprintf('TrendInflationData-%s-%s', Ycode{n}, datalabel), wrap);
end

%% simulate predictive density at end of sample
NfcstDraws    = 1e2; % draws per node

fcstYdraws = NaN(N, NfcstDraws, NfcstHorizons, Nmcmc); % note: will be permuted below
C          = [ones(N,1) eye(N,N*p)];

for mm = 1 : Nmcmc


    PAI          = drawPAI(:,:,mm);
    Y0           = transpose(drawY0(:,:,mm));
    StateBar     = transpose(drawStateBar(end,:,mm));
    StateGap     = transpose(drawStateGap(end,:,mm));
    sqrtOMEGAgap = drawSQRTOMEGAgap(:,:,mm);
    x0           = repmat(drawXjumpoff(:,mm), [1 NfcstDraws]);

    logSV0       = drawlogSVbar(:,end,mm);
    varlogSV     = drawlogvarSVbar(:,mm);

    B                = zeros(Nbar+Ngap*p,Nbar+Ngap);
    B(ndxbar,ndxbar) = eye(Nbar);
    B(ndxgap,ndxgap) = sqrtOMEGAgap;

    SVshocks       = randn(rndStream, Nbar, NfcstDraws, NfcstHorizons);
    SVpath         = logSV0 + varlogSV .* cumsum(SVshocks);
    SVpath         = exp(.5 * SVpath);

    shocks             = randn(rndStream, Nbar + Ngap, NfcstDraws, NfcstHorizons);
    shocks(ndxbar,:,:) = shocks(ndxbar,:,:) .* SVpath;

    gapkompanion                    = zeros(N * p);
    gapkompanion(N+1:end,1:N*(p-1)) = eye(N*(p-1));
    gapkompanion(1:N,:)             = PAI;
    kompanion                       = blkdiag(1, gapkompanion);

    for hh = 1 : NfcstHorizons
        x                     = kompanion * x0 + B * shocks(:,:,hh);
        fcstYdraws(:,:,hh,mm) = C * x + Y0;
        x0                    = x;
    end
end

fcstYdraws = permute(fcstYdraws, [1 3 2 4]);
fcstYdraws = reshape(fcstYdraws, [N, NfcstHorizons, NfcstDraws * Nmcmc]);

%% plot predictive density
fcstYhat   = mean(fcstYdraws, 3);
fcstYtail  = prctile(fcstYdraws, normcdf([-1 1]) * 100, 3);

fcstDates  = dateshift(dates(end), 'start', 'month', 1:NfcstHorizons);
for n = 1 : N
    thisfig = figure;
    hold on
    plot(fcstDates,  fcstYhat(n,:), 'k-', 'linewidth', 2);
    plot(fcstDates,  squeeze(fcstYtail(n,:,:)), 'k-', 'linewidth', 1);
    yline(medTrends(end,n), '-.', 'color', Colors4Plots('darkred'), 'linewidth', 3)
    % yline(tailTrends(end,n,1), '-.', 'color', Colors4Plots('darkred'), 'linewidth', 1)
    % yline(tailTrends(end,n,2), ':', 'color', Colors4Plots('darkred'), 'linewidth', 1)
    xlim(fcstDates([1 end]))
    ylim0
    grid on
    title(Ycode{n})
    wrapthisfigure(thisfig, sprintf('predictivedensity-%s-%s', Ycode{n}, datalabel), wrap);
end

%% finish
dockAllFigures
finishwrap
if ~isempty(wrap) && (~isdesktop || ispc)
    close all
    thisArch = computer('arch');
    thisVer  = ver;
    if ismac
        [~, thisSys]   = system('sysctl -a | grep machdep.cpu ', '-echo');
        [~, thisBrand] = system('sysctl -a | grep machdep.cpu | grep brand_string ', '-echo');
        if contains(thisBrand, 'M1 Pro')
            thisBrand = 'AppleSilicon';
        else
            thisBrand = 'MacOSIntel';
        end
    elseif isunix
        [~, thisSys] = system('cat /proc/cpuinfo ', '-echo');
        thisBrand = 'IntelUbuntu';
    else % ispc
        thisSys   = 'Intel(R) Xeon(R) Gold 6320 CPU @ 2.1 GHz';
        thisBrand = 'WindowsXeon';
    end
    varlist = {'thisArch', 'thisVer', 'thisSys', 'thisBrand', ...
        'MCMCtime', 'execLabel', 'doSingleThread', 'usedThreads', 'availableThreads', ...
        'datalabel', 'Nmcmc', 'burnin', 'dates', ...
        'Yinflation', 'N', 'Ycode', 'np', ...
        'fcstYhat', 'fcstYtail', 'fcstDates', ...
        'med*', '*tail*'};
    save(fullfile(wrap.dir, strcat('slim-', titlename)), varlist{:}, '-v7.3')
    % save(fullfile(wrap.dir, titlename), '-v7.3')
end
finishscript