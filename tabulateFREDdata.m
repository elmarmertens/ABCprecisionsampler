%% load FRED-based data sets and tabulate availability

clear
clc
close all


datalabel = 'INFTRMSRV';

FREDtable = readtable(sprintf('%s.csv', datalabel));

dates     = datetime(FREDtable.dates, 'ConvertFrom', 'datenum');
VarNames  = FREDtable.Properties.VariableNames(2:end);
data      = table2array(FREDtable(:,2:end));
N = length(VarNames);

initwrap

%% collect data
tabheaders = {'Variable', 'Frequency', 'Avail. since', 'Number of obs.'};

tabfreq                             = cell(N,1);
tabfreq(:)                          = {'monthly'};
tabfreq(contains(VarNames, 'SPF'))  = {'quarterly'};
tabfreq(contains(VarNames, 'GDPD')) = {'quarterly'};

tabobs1 = cell(N,1);
for n = 1 : N
    t1         = find(~isnan(data(:,n)), 1);
    tabobs1{n} = dates(t1);
end

tabnobs = sum(~isnan(data), 1);

%% generate pretty labels
prettyVarNames = cell(size(VarNames));
for n = 1 : N
    switch VarNames{n}
        case 'PCE'
            prettyVarNames{n} = 'headline PCE';
        case 'PCEcore'
            prettyVarNames{n} = 'core PCE';
        case 'CPI'
            prettyVarNames{n} = 'headline PCE';
        case 'GDPD'
            prettyVarNames{n} = 'GDP deflator';
        case 'PCEtrim'
            prettyVarNames{n} = 'trimmed PCE (Dallas Fed)';
        case 'CPItrim'
            prettyVarNames{n} = 'trimmed CPI (Cleveland Fed)';
        case 'CPImedian'
            prettyVarNames{n} = 'median CPI (Cleveland Fed)';
        case 'CPIsticky'
            prettyVarNames{n} = 'sticky CPI (headline, Atlanta Fed)';
        case 'CPIcoresticky'
            prettyVarNames{n} = 'sticky CPI (core, Atlanta Fed)';
        case 'SPFcpi10Y'
            prettyVarNames{n} = '10-year expected CPI (SPF)';
        case 'SPFcpi1Y'
            prettyVarNames{n} = 'SPF: CPI change over next four quarters';
        case 'SPFpce1Y'
            prettyVarNames{n} = 'SPF: PCE change over next four quarters';
        case 'SPFgdpd1Y'
            prettyVarNames{n} = 'SPF: GDP deflator over next four quarters';
        otherwise
            prettyVarNames{n} = VarNames{n};
    end
end
%% create table
tabname = sprintf('datalist%s.tex', datalabel);
tabdir  = wrap.dir;
latexwrapper(wrap, 'add', 'tab', tabname)

fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
Ncols = length(tabheaders);
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('c', 1, Ncols));
fprintf(fid, '\\toprule\n');

% column headers
fprintf(fid, '%s', tabheaders{1});
fprintf(fid, ' & %s', tabheaders{2:end});
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
% loop over rows
for n = 1 : N
    fprintf(fid, '%s ', prettyVarNames{n});
    fprintf(fid, '& %s ', tabfreq{n});
    fprintf(fid, '& %s ', datetime(tabobs1{n}, 'Format', 'yyyy MMM'));
    fprintf(fid, '& %d ', tabnobs(n));
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');

fprintf(fid, 'Notes: ');
fprintf(fid, 'All data has been obtained from the St. Louis Fed''s FRED database per May 9 2023.\n');
fprintf(fid, '(For realized inflation measures, seasonally adjusted series have been obtained.)\n');
fprintf(fid, 'Unless otherwise noted, headline, core, sticky, trimmed and median inflation are measured as monthly changes. GDP deflator inflation as quarterly change.\n');
fprintf(fid, 'Trimmed PCE inflation and sticky core inflation measures reflect 12-month changes.\n');
fprintf(fid, 'Each measure is expressed as a continuously compounded, annualized rate of inflation.\n');
fprintf(fid, '\n');

fclose(fid);
type(fullfile(tabdir, tabname))

%% finish
finishwrap