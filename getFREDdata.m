%% construct FRED-based data sets
% uses input files obtained from FRED and FRB-PHIL
%#ok<*DATNM>

clear
clc
close all

samStart = datenum(1960,1,1);
samEnd   = datenum(2023,4,1);

dates    = genrMdates(1960,year(samEnd),1);
dates    = dates(dates <= samEnd);
T        = length(dates);

DATALABEL = {'INF', 'INFTRM', 'INFTRMSRV'};

freddir = 'INFTRMINT'; %note: need to manually remove "." entries for missing obs of GS20 from FRED file

%% process datalabel
for d = 1 : length(DATALABEL)

    datalabel = DATALABEL{d};


    % mndx1, and qndx1 refer to columns in data output
    % mndx2 refers to column in FRED source file

    switch datalabel
        case 'INF'
            Ylabel   = {'PCE', 'PCEcore', 'CPI', 'GDPD'};
            mndx1    = [1 2 3];
            mndx2    = [12 13 3];
            qndx1    = 4;
            useSPF   = false;
        case 'INFTRMSRV'
            Ylabel  = {'PCE', 'PCEcore', 'CPI', 'GDPD', ...
                'PCEtrim', 'CPItrim', 'CPImedian', ...
                'CPIsticky', 'CPIcoresticky', ...
                'SPFcpi10Y', 'SPFcpi1Y', 'SPFpce1Y', 'SPFgdpd1Y'}; 
            mndx1   = [1:3, 5:9];
            mndx2   = [12 13 3 14 16 11 15 2];
            qndx1   = 4;
            spfndx1 = 10:13;
            useSPF = true;
        otherwise
            error('datalabel <<%s>> not known', datalabel)
    end

    Ny = length(Ylabel);

    %% load FRED data

    % qd = importdata(fullfile(freddir, sprintf('%s_Quarterly.txt', freddir)));
    % md = importdata(fullfile(freddir, sprintf('%s_Monthly.txt', freddir)));
    qd     = readtable(fullfile(freddir, sprintf('%s_Quarterly.txt', freddir)));
    qData  = table2array(qd(:,2:end));
    qDates = datenum(qd.DATE);

    md     = readtable(fullfile(freddir, sprintf('%s_Monthly.txt', freddir)));
    mData  = table2array(md(:,2:end));
    mDates = datenum(md.DATE);

    %% transform all series into annualized log-changes
    qData(2:end,:) = diff(log(qData)) * 400;
    qData(1,:)     = NaN;


    % COREFLEXCPIM159SFRBATL	CORESTICKM159SFRBATL	CPIAUCSL	CPILFESL	GS10	GS2	GS20	GS30	GS5	GS7	MEDCPIM158SFRBCLE	PCEPI	PCEPILFE	PCETRIM12M159SFRBDAL	STICKCPIM157SFRBATL	TRMMEANCPIM158SFRBCLE
    levelNdx = [3 4 12 13];
    ppAPRNdx = [1 2 11 14 16 5:10];
    ppMMndx  = 15; % sticky CPI

    mData(2:end,levelNdx) = diff(log(mData(:,levelNdx))) * 1200;
    mData(1,levelNdx)     = NaN;

    mData(:,ppAPRNdx)     = log(1 + mData(:,ppAPRNdx) / 100) * 100;
    mData(:,ppMMndx)      = log(1 + mData(:,ppMMndx) / 100) * 1200;




    %% prep destintion data and copy FRED data

    data = NaN(T,Ny);

    m2dates = ismember(mDates, dates);
    data(ismember(dates, mDates), mndx1) = mData(m2dates,mndx2);

    q2dates = ismember(qDates, dates);
    if ~isempty(qndx1)
        data(ismember(dates, qDates), qndx1) = qData(q2dates,:);
    end
    % push quarterly data two months later
    data(3:end,qndx1) = data(1:end-2,qndx1);
    data(1:2,qndx1)   = NaN;

    %% load SPF data
    if useSPF

        Nspf     = 1;

        % first, read in CPI10Y
        CPI10Y   = importdata('Median_CPI10_Level.xlsx');
        spfY     = CPI10Y.data(:,1);
        spfQ     = CPI10Y.data(:,2);
        spfdates = datenum(spfY, (spfQ - 1) * 3 + 2, 1); % assign to mid-of-quarter month

        SPFdata  = NaN(length(spfdates), Nspf);
        SPFdata(:,1) = CPI10Y.data(:,3);

        % read PCE, CPI etc.
        SPFlabels = {'CPI', 'PCE', 'PGDP'};
        for s = 1 : length(SPFlabels)
            spfimport = importdata(sprintf('Median_%s_Level.xlsx', SPFlabels{s}));
            if any(spfY ~= CPI10Y.data(:,1)) || any(spfQ  ~= CPI10Y.data(:,2))
                error('date mismatch')
            end
            switch SPFlabels{s}
                case 'PGDP'
                    thisdata = (log(spfimport.data(:,8)) - log(spfimport.data(:,4))) * 100;
                    thisdata = mean(thisdata,2);
                otherwise
                    thisdata = log(1 + spfimport.data(:,5:8) / 100) * 100;
                    thisdata = mean(thisdata,2);
            end
            SPFdata(:,1+s) = thisdata;

        end

        % patch SPFdata into data
        [~, ndx1, ndx2]   = intersect(dates, spfdates);
        data(ndx1,spfndx1) = SPFdata(ndx2,:);

    end


    %% write output files
    ynan = isnan(data);
    data(ynan) = 0;
     
    % mat2fortran(sprintf('%s.dates.txt', datalabel), dates);
    % logical2fortran(sprintf('%s.yNaN.txt', datalabel), ynan);
    % mat2fortran(sprintf('%s.yData.txt', datalabel), data);
    % 
    % filename = sprintf('%s.settings.txt', datalabel);
    % fid = fopen(filename, 'wt');
    % fprintf(fid, 'Ny = %d\n', size(data,2));
    % fprintf(fid, 'T  = %d\n', size(data,1));
    % fprintf(fid, 'YLABEL:\n');
    % for n = 1 : Ny
    %     fprintf(fid, '%s\n', Ylabel{n});
    % end
    % fclose(fid);
    % display(filename);
    % type(filename)
    % hrulefill
    % 
    % % save also as matfile
    % varlist = {'dates', 'ynan', 'data', 'datalabel', 'Ylabel'};
    % save(datalabel, varlist{:}, '-v7.3')

    % CSV TABLE
    data(ynan) = NaN;
    FREDtable = array2table([dates data], 'VariableNames', ['dates', Ylabel]);
    writetable(FREDtable, sprintf('%s.csv', datalabel))

    %% plot
    figure
    plot(dates, data)
    xtickdates(dates)

end % for d
