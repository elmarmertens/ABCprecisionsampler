%% collect matfiles and report MCMCtoc times

clear
clc
close all

matdir    = 'matfilesJEDCR1';
datalabel = 'INFTRMSRV';
Nmcmc     = 1e4;
initwrap

%% collect matfiles
mates = dir(fullfile(matdir, sprintf('slim-commontrendinflation-%s-Ndraws%s-*multithreaded.mat', datalabel,  erase(sprintf('%1.0e', Nmcmc), '+'))));
% mates = dir(fullfile(matdir, 'commontrendinflationRbar-*.mat'));

matefiles = cell(length(mates), 1);

for m = 1 : length(mates)
    matefiles{m} = matfile(fullfile(mates(m).folder, mates(m).name));
end

datalabels  = cellfun(@(c) c.datalabel, matefiles, 'UniformOutput',false);
execLabels  = cellfun(@(c) c.execLabel, matefiles, 'UniformOutput',false);
MCMCseconds = cellfun(@(c) c.MCMCtime, matefiles, 'UniformOutput',true);
MCMCminutes = MCMCseconds / 60;
MCMChours   = MCMCminutes / 60;

MCMCrem     = mod(MCMCminutes, 60);

%% construct pretty label
prettyLabels = cell(size(execLabels));
for nn = 1 : length(execLabels)
    if contains(execLabels(nn), 'DK')
        prettyLabels{nn} = 'DK';
    else
        prettyLabels{nn} = 'PS';
    end
    % if contains(execLabels(nn), 'multi')
    %     prettyLabels{nn} = strcat(prettyLabels{nn}, ' (multi-threaded)');
    % else
    %     prettyLabels{nn} = strcat(prettyLabels{nn}, ' (single thread)');
    % end
end

%% tabulate
timetable = table(datalabels, execLabels, prettyLabels, MCMCseconds, MCMCminutes, MCMChours, MCMCrem);
timetable = sortrows(timetable, [1 4], 'descend');
display(timetable)

%% barplot

fontsize = 14;
thisfig = figure;
bar(1:length(mates), timetable.MCMChours)
ax = gca;
set(ax, 'fontsize', fontsize)
xticklabels(timetable.prettyLabels)
ylim([0 15])
valLabels = arrayfun(@(h,m) sprintf('%3.0fh %dmin', h, m), floor(timetable.MCMChours), round(timetable.MCMCrem), 'uniformoutput', false);
text(1:length(mates),timetable.MCMChours, valLabels,'vert','bottom','horiz','center', 'Fontsize', fontsize); 
box off
wrapthisfigure(thisfig, 'commontrendMCMCtimesMultithreaded', wrap)

%% finish
finishwrap
