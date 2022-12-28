%% tabulate performance times for VAR's with missing value

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


%% load results

for BRANDS = {'MacOSIntel', 'WindowsXeon'}

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

        for missValShare = [.01 .05 .1 .2 .3 .5]

            if strcmp(thisBrand, 'IntelUbuntu') && missValShare > .01
                continue % skip other cases, since matfiles for IntelUbuntu are missing
            end

            matfiledir = 'matfiles2';
            matname = sprintf('VARmissingvaluesPStimes%sThreads%dof%dmissValShare%02d', thisBrand, usedThreads, availableThreads, floor(missValShare * 100));
            mat     = matfile(fullfile(matfiledir, matname));

            if mat.missValShare ~= missValShare
                error('missValShare mismatch in matfile')
            end

            thisVer    = mat.thisVer; % MatFile object does not support indexing into struct arrays
            ndx        = arrayfun(@(x) strcmp(x.Name, 'MATLAB'), thisVer);
            thisMatlab = sprintf('%s %s', thisVer(ndx).Name, thisVer(ndx).Release);


            gridNy = mat.gridNy;
            gridT  = mat.gridT;
            gridP  = mat.gridP;

            % [PStimes, PS0times, DKtimes]  = deal(NaN(length(gridP), length(gridNy), length(gridT)));
            PStimes      = mat.PStimes;
            PS0times     = mat.PS0times;
            DKtimes      = mat.DKtimes;
            DKMisstimes  = mat.DKMisstimes;

            if length(gridT) ~= 2
                error('expecting gridT to have only 2 elements');
            end


            % PStimes  ./ DKtimes
            % PS0times ./ DKtimes


            %% create latex table

            panellabel   = {'A', 'B', 'C', 'D'};
            panelcaption = {'PS (excl. QR) relative to DK (in pp)', 'PS (incl. QR) relative to DK (in pp)', 'DK-missing-variables relative to DK (in pp)', 'DK (in seconds)'};
            Npanel       = length(panellabel);

            tabname    = sprintf('VARmissingvaluesPStimes%sThreads%dof%dmissValShare%02d.tex', thisBrand, usedThreads, availableThreads, floor(missValShare * 100));
            %  tabcaption = sprintf('Execution times for sampling missing values of a VAR');
            tabcaption = sprintf('Missing-value VAR on %s (%s) with %d\\%% of observations missing', thisComputerShort, threadLabel, floor(missValShare * 100));
            tabdir     = wrap.dir;

            latexwrapper(wrap, 'add', 'tab', tabname, tabcaption)

            % collect table data
            tableData = NaN(length(gridP), length(gridNy) * 2, Npanel);
            iterPanel = 1;
            thisPanel = PStimes ./ DKtimes * 100;
            tableData(:,:,iterPanel) = reshape(thisPanel, [length(gridP), length(gridNy) * 2]);
            iterPanel = 2;
            thisPanel = PS0times ./ DKtimes * 100;
            tableData(:,:,iterPanel) = reshape(thisPanel, [length(gridP), length(gridNy) * 2]);
            iterPanel = 3;
            thisPanel = DKMisstimes ./ DKtimes * 100;
            tableData(:,:,iterPanel) = reshape(thisPanel, [length(gridP), length(gridNy) * 2]);
            iterPanel = 4;
            thisPanel = DKtimes;
            tableData(:,:,iterPanel) = reshape(thisPanel, [length(gridP), length(gridNy) * 2]);

            % tabulate
            fid = fopen(fullfile(tabdir, tabname), 'wt');
            fprintf(fid, '\\begin{center}\n');
            Ncols = 1 + 2 * length(gridNy);
            fprintf(fid, '\\begin{tabular}{r%s}\n', repmat('.1', 1, 2 * length(gridNy)));
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
                        fprintf(fid, '& %8.0f ', tableData(iterP,:,iterPanel));
                    else
                        fprintf(fid, '& \\multicolumn{1}{r}{%8.2f} ', tableData(iterP,:,iterPanel));
                    end
                    fprintf(fid, '\\\\\n');
                end
            end

            fprintf(fid, '\\bottomrule\n');
            fprintf(fid, '\\end{tabular}\n');
            fprintf(fid, '\\end{center}\n');
            fprintf(fid, '\n');

            fprintf(fid, 'Notes: ');
            fprintf(fid, 'Based on simulated data with %d\\%% of all observations missing.\n', floor(missValShare * 100));
            fprintf(fid, 'Panels~A and~B report the execution time of a typical call to the precision-based sampler (PS) for different choices of lag length ($p$),  number of VAR variables ($N_y$) and observations ($T$) in percentage points of the execution time of the Durbin-Koopmann''s disturbance smoothing sampler (DK) whose execution time (in seconds) is reported in Panel~D.\n');
            fprintf(fid, 'Execution times for the  precision-based sampler reported in Panel A reflect the use of prepared one-off computations (incl. the QR decomposition of measurement loadings $\\boldsymbol{C}$) outside the measured times. These time are relevant for MCMC applications where $\\boldsymbol{C}$ does not change between sampling steps.\n');
            fprintf(fid, 'Panel~B considers calls to the precision-based sampler that encompass all computations.\n');
            fprintf(fid, 'Panel~C considers a version of the DK sampler that is specialized to the missing-value case (where the measurement loadings merely pick specific rows of the state vector).\n');
            fprintf(fid, 'All times were measured in %s with the \\texttt{timeit} function on an %s ', thisMatlab, thisComputer);
            if usedThreads == 1
                fprintf(fid, ' in single-threaded mode.\n');
            else
                fprintf(fid, ' in multi-threaded mode (for matrix operations) using %d threads.\n', usedThreads);
            end
            fprintf(fid, '\n');

            fclose(fid);
            type(fullfile(tabdir, tabname))

        end % missValShare
    end % usedThreads
end % thisBrand

%% finish
finishwrap
finishscript