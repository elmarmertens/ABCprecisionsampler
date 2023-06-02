%% plot trend inflation estimates


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


resultsdir = 'matfilesJEDCR1';

%% select model

datalabel = 'INFTRMSRV';
Nmcmc     = 1e4;
Nbar      = 1;

execLabel      = 'PS';
doSingleThread = false;

titlename = sprintf('commontrendinflation-%s-Ndraws%s', datalabel, erase(sprintf('%1.0e', Nmcmc), '+'));
if doSingleThread
    execLabel = strcat(execLabel, 'singlethreaded');
else
    execLabel = strcat(execLabel, 'multithreaded');
end
titlename = strcat(titlename, '-', execLabel);

%% load results
varlist = {'thisArch', 'thisVer', 'thisSys', 'thisBrand', ...
    'MCMCtime', 'execLabel', 'doSingleThread', 'usedThreads', 'availableThreads', ...
    'datalabel', 'Nmcmc', 'burnin', 'dates', ...
    'Yinflation', 'N', 'Ycode', 'np', ...
    'fcstYhat', 'fcstYtail', 'fcstDates', ...
    'med*', '*tail*'};
load(fullfile(resultsdir, strcat('slim-', titlename)), varlist{:})


%% start output
initwrap

fprintf('MCMC took %8.2f (with %s)\n', MCMCtime, execLabel)


%% trend variance
for n = 1 : Nbar
    thisfig = figure;
    hold on
    plot(dates, medSV(n,:), 'k-', 'LineWidth', 2);
    plot(dates, squeeze(tailSV(n,:,:)), 'k-', 'LineWidth', 1);
    ylim0
    wrapthisfigure(thisfig, sprintf('SVbar%d-%s',  n, datalabel), wrap);
end


%% report Trends
fontsize = 16;
YLIM = [-2 16];
for n = 1 : N
    thisfig = figure;
    hold on
    ax = gca;
    set(ax, 'fontsize', fontsize)
    plotCI(medY(:,n), squeeze(tailY(:,n,:)), dates, [], ':', 'color', [.25 .25 .25], 'linewidth', 2);
    h = plot(dates,  medTrends(:,n), '-', 'linewidth', 2, 'color', Colors4Plots('darkred'));
    plot(dates,  squeeze(tailTrends(:,n,:)), '-', 'linewidth', 1, 'color', Colors4Plots('darkred'));
    if iscompact(Yinflation(:,n))
        hData = plot(dates,Yinflation(:,n), ':', 'color', 'black', 'linewidth', 1);
    else
        hData = plot(dates,Yinflation(:,n), 'd', 'color', 'black', 'linewidth', 1);
    end
    % hDataMA = plot(dates, meanK(Yinflation(:,n), np, 1), '-', 'color', 'black', 'linewidth', 1);
    xtickdates(dates)
    ylim0
    % ylim([-2 max(ylim)])
    grid on
    % title(Ycode{n})
    wrapthisfigure(thisfig, sprintf('TrendInflationData-%s-%s', Ycode{n}, datalabel), wrap);
    legend([h hData], 'Trend', 'Data') % , sprintf('Data MA(%d)', np))
    wrapthisfigure(thisfig, sprintf('TrendInflationData-%s-%s-WITHLEGEND', Ycode{n}, datalabel), wrap);
end

%% compare against INF trend
mat = matfile(fullfile(resultsdir, 'slim-commontrendinflation-INF-Ndraws1e03-PSmultithreaded.mat'));
for n = 1 : 4
    thisfig = figure;
    hold on
    ax = gca;
    set(ax, 'fontsize', fontsize)
    % plotCI(medY(:,n), squeeze(tailY(:,n,:)), dates, [], ':', 'color', [.25 .25 .25], 'linewidth', 2);
    h2 = plot(dates,  mat.medTrends(:,n), '-.', 'linewidth', 2, 'color', Colors4Plots('green'));
    h = plot(dates,  medTrends(:,n), '-', 'linewidth', 2, 'color', Colors4Plots('darkred'));
    plot(dates,  squeeze(tailTrends(:,n,:)), '-', 'linewidth', 1, 'color', Colors4Plots('darkred'));
    % if iscompact(Yinflation(:,n))
    %     hData = plot(dates,Yinflation(:,n), ':', 'color', 'black', 'linewidth', 1);
    % else
    %     hData = plot(dates,Yinflation(:,n), 'd', 'color', 'black', 'linewidth', 1);
    % end
    hDataMA = plot(dates, meanK(Yinflation(:,n), np, 1), '-', 'color', 'black', 'linewidth', 1);
    xtickdates(dates)
    ylim0
    % ylim([-2 max(ylim)])
    grid on
    % title(Ycode{n})
    wrapthisfigure(thisfig, sprintf('TrendInflation2Data-%s-%s', Ycode{n}, datalabel), wrap);
    legend([h h2 hData], 'Trend', 'Trend (INF)', 'Data') % , sprintf('Data MA(%d)', np))
    wrapthisfigure(thisfig, sprintf('TrendInflation2Data-%s-%s-WITHLEGEND', Ycode{n}, datalabel), wrap);
end


%% plot predictive density
% fcstDates  = dateshift(dates(end), 'start', 'month', 1:NfcstHorizons);
% for n = 1 : N
%     thisfig = figure;
%     hold on
%     plot(fcstDates,  fcstYhat(n,:), 'k-', 'linewidth', 2);
%     plot(fcstDates,  squeeze(fcstYtail(n,:,:)), 'k-', 'linewidth', 1);
%     yline(medTrends(end,n), '-.', 'color', Colors4Plots('darkred'), 'linewidth', 3)
%     % yline(tailTrends(end,n,1), '-.', 'color', Colors4Plots('darkred'), 'linewidth', 1)
%     % yline(tailTrends(end,n,2), ':', 'color', Colors4Plots('darkred'), 'linewidth', 1)
%     xlim(fcstDates([1 end]))
%     ylim0
%     grid on
%     title(Ycode{n})
%     wrapthisfigure(thisfig, sprintf('predictivedensity-%s-%s', Ycode{n}, datalabel), wrap);
% end


%% finish
dockAllFigures
finishwrap
finishscript