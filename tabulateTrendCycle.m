%% tabulate performance times for Trend-cum-VAR model

% coded by Elmar Mertens em@elmarmertens.com


%% clean workspace
clear
clc
close all

%% load toolboxes
path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%#ok<*UNRCH>
%#ok<*NOPTS>

%% prepare latex wrapper

initwrap

matfiledir = 'matfilesJEDCR1';

VARLABEL   = {'COMMONTRENDCYCLE', 'TRENDCYCLE' };
VARCAPTION = {'Common-trend model', 'Trend-cycle-VAR model'};


%% load results

for vv = 1  : length(VARLABEL)

    VARlabel   = VARLABEL{vv};
    VARcaption = VARCAPTION{vv};

    for BRANDS = {'WindowsXeon', 'AppleSilicon', 'MacOSIntel'}

        thisBrand = BRANDS{:};
        switch thisBrand
            case 'WindowsXeon'
                thisComputerShort = 'Intel Xeon w/Windows';
                thisComputer = 'Intel(R) Xeon(R) Gold 6320 CPU @ 2.1 GHz (Windows)';
                availableThreads = 32;
            case 'MacOSIntel'
                thisComputerShort = 'Intel i7 w/macOS';
                thisComputer = 'Intel(R) Core(TM) i7-5775R CPU @ 3.30GHz (macOS)';
                availableThreads = 4;
            case 'AppleSilicon'
                thisComputerShort = 'Apple Silicon';
                thisComputer = 'Apple M1 Pro (macOS, Rosetta 2)';
                availableThreads = 8;
            case 'IntelUbuntu'
                thisComputerShort = 'Intel Xeon w/Ubuntu';
                thisComputer = 'Intel(R) Xeon(R) CPU E5-1680 v4 @ 3.40GHz';
                availableThreads = 8;
            otherwise
                error('brand <<%s>> not known.', thisBrand)
        end

        for usedThreads = [1 availableThreads]

            if usedThreads == 1
                threadLabel = sprintf('%d thread', usedThreads);
            else
                threadLabel = sprintf('%d threads', usedThreads);
            end

            if strcmp(thisBrand, 'IntelUbuntu') && missValShare > .01
                continue % skip other cases, since matfiles for IntelUbuntu are missing
            end

            matname = sprintf('%sPStimes%sThreads%dof%d', VARlabel, thisBrand, usedThreads, availableThreads);
            mat     = matfile(fullfile(matfiledir, matname));

            thisVer    = mat.thisVer; % MatFile object does not support indexing into struct arrays
            ndx        = arrayfun(@(x) strcmp(x.Name, 'MATLAB'), thisVer);
            thisMatlab = sprintf('%s %s', thisVer(ndx).Name, thisVer(ndx).Release);


            gridNy = mat.gridNy;
            gridT  = mat.gridT;
            gridP  = mat.gridP;

            % [PStimes, PS0times, DKtimes]  = deal(NaN(length(gridP), length(gridNy), length(gridT)));
            PSyxtimes  = mat.PSyxtimes;
            TCPStimes  = mat.TCPStimes;
            TCPS0times = mat.TCPS0times;
            PStimes    = mat.PStimes;
            PS0times   = mat.PS0times;
            DKtimes    = mat.DKtimes;
            PSnoisetimes  = mat.PSnoisetimes;

            if length(gridT) ~= 2
                error('expecting gridT to have only 2 elements');
            end


            % PStimes  ./ DKtimes
            % PS0times ./ DKtimes


            %% create latex table

            panellabel   = {'A', 'B', 'C', 'D', 'E'};
            panelcaption = {'ABC-PS (excl. QR) as percentage of DK', 'ABC-PS (incl. QR) as percentage of DK', ...
                'Trend-Cycle PS as percentage of DK', ...
                'standard PS w/noise as percentage of DK', ...
                'DK in seconds'};
            Npanel       = length(panellabel);

            tabname    = sprintf('%sPStimes%sThreads%dof%d.tex', VARlabel, thisBrand, usedThreads, availableThreads);
            %  tabcaption = sprintf('Execution times for sampling missing values of a VAR');
            tabcaption = sprintf('%s on %s (%s)', VARcaption, thisComputerShort, threadLabel);
            tabdir     = wrap.dir;

            latexwrapper(wrap, 'add', 'tab', tabname, tabcaption)

            % collect table data
            tableData = NaN(length(gridP), length(gridNy) * 2, Npanel);
            iterPanel = 0;
            iterPanel = iterPanel + 1;
            thisPanel = PStimes ./ DKtimes * 100;
            tableData(:,:,iterPanel) = reshape(thisPanel, [length(gridP), length(gridNy) * 2]);
            iterPanel = iterPanel + 1;
            thisPanel = PS0times ./ DKtimes * 100;
            tableData(:,:,iterPanel) = reshape(thisPanel, [length(gridP), length(gridNy) * 2]);
            iterPanel = iterPanel + 1;
            thisPanel = TCPStimes ./ DKtimes * 100;
            tableData(:,:,iterPanel) = reshape(thisPanel, [length(gridP), length(gridNy) * 2]);
            iterPanel = iterPanel + 1;
            thisPanel = PSnoisetimes ./ DKtimes * 100;
            tableData(:,:,iterPanel) = reshape(thisPanel, [length(gridP), length(gridNy) * 2]);
            iterPanel = iterPanel + 1;
            thisPanel = DKtimes;
            tableData(:,:,iterPanel) = reshape(thisPanel, [length(gridP), length(gridNy) * 2]);

            % tabulate
            fid = fopen(fullfile(tabdir, tabname), 'wt');
            fprintf(fid, '\\begin{center}\n');
            Ncols = 1 + 2 * length(gridNy);
            fprintf(fid, '\\begin{tabular}{r%s}\n', repmat('.2', 1, 2 * length(gridNy)));
            fprintf(fid, '\\toprule\n');

            % column headers
            fprintf(fid, ' & ');
            fprintf(fid, '\\multicolumn{%d}{c}{$N_y$} ', 2 * length(gridNy));
            fprintf(fid, '\\\\\n');
            fprintf(fid, '\\cmidrule(lr){%d-%d}\n',2,Ncols);


            iterT = 1;
            fprintf(fid, '& \\multicolumn{%d}{c}{$T = %d$} ', length(gridNy), gridT(iterT));
            iterT = 2;
            fprintf(fid, '& \\multicolumn{%d}{c}{$T = %d$} ', length(gridNy), gridT(iterT));
            fprintf(fid, '\\\\\n');
            fprintf(fid, '\\cmidrule(lr){%d-%d}\n',2,1 + length(gridNy));
            fprintf(fid, '\\cmidrule(lr){%d-%d}\n',2 + length(gridNy), Ncols);

            fprintf(fid, '$p$ ');
            fprintf(fid, '& \\multicolumn{1}{c}{\\quad$%d$\\quad} ', gridNy);
            fprintf(fid, '& \\multicolumn{1}{c}{\\quad$%d$\\quad} ', gridNy);
            fprintf(fid, '\\\\\n');

            for iterPanel = 1 : Npanel
                fprintf(fid, '\\midrule\n');
                fprintf(fid, '\\multicolumn{%d}{c}{\\textbf{PANEL %s: %s}} ', Ncols, panellabel{iterPanel}, panelcaption{iterPanel});
                fprintf(fid, '\\\\\n');
                fprintf(fid, '\\midrule\n');
                for iterP = 1 : length(gridP)
                    fprintf(fid, '$%d$ ', gridP(iterP));
                    if iterPanel < Npanel
                        fprintf(fid, '& \\multicolumn{1}{r}{%4.0f} ', tableData(iterP,:,iterPanel));
                    else
                        fprintf(fid, '& %8.2f ', tableData(iterP,:,iterPanel));
                    end
                    fprintf(fid, '\\\\\n');
                end
            end

            fprintf(fid, '\\bottomrule\n');
            fprintf(fid, '\\end{tabular}\n');
            fprintf(fid, '\\end{center}\n');
            fprintf(fid, '\n');

            fprintf(fid, 'Notes: ');
            fprintf(fid, 'Based on simulated data.\n');
            fprintf(fid, 'Panels~A through~D report execution times of precision-based samplers (PS) as percentage of the execution time of the Durbin-Koopmann''s disturbance smoothing sampler (DK).\n');
            fprintf(fid, 'Execution times (in seconds) of the DK sampler are reported in Panel~E.\n');
            fprintf(fid, 'Panels~A and~B report execution times for the generic ABC precision-based sampler, with Panel~B reflecting the use of prepared one-off computations.\n');
            fprintf(fid, 'Panel~C provides results for a PS specicalized to the trend-cycle case (and prepared one-off computations).\n');
            fprintf(fid, 'Panel~D reflects calls to a standard precision-based sampler, when called with a minimal noise term added to the measurement equation.\n');
            fprintf(fid, 'All times were measured in %s with the \\texttt{timeit} function on an %s ', thisMatlab, thisComputer);
            if usedThreads == 1
                fprintf(fid, ' in single-threaded mode.\n');
            else
                fprintf(fid, ' in multi-threaded mode (for matrix operations) using %d threads.\n', usedThreads);
            end
            fprintf(fid, '\n');

            fclose(fid);
            type(fullfile(tabdir, tabname))

        end % usedThreads
    end % thisBrand
end % VARlabel

%% finish
finishwrap
finishscript