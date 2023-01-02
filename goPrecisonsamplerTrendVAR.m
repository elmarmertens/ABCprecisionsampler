%% Apply precision-based ABC sampler to Common-Trend-cum-VAR(p) models (and simulated data)


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

rng(061222); %  fix random seed

quicky = false;

%#ok<*UNRCH>
%#ok<*NOPTS>

rng(011222)

for doSingleThread = [true false]

    maxNumCompThreads('automatic');
    if doSingleThread
        % enforce single threaded compuations
        usedThreads = 1;
        availableThreads = maxNumCompThreads(usedThreads);
        fprintf('Using 1 of %d available threads.\n', availableThreads)
    else
        usedThreads =  maxNumCompThreads('automatic');
        availableThreads = usedThreads;
        fprintf('Using all %d available threads.\n', usedThreads)
    end


    %% define parameter grid

    gridNy = 5 : 5 : 25;
    gridT  = [200 800];
    gridP  = [4 8 12 24];

    [PStimes, PS0times, PSyxtimes, PS0yxtimes, DKtimes]  = deal(NaN(length(gridP), length(gridNy), length(gridT)));

    for iterT = 1 : length(gridT)
        T = gridT(iterT)
        for iterNy = 1 : length(gridNy)
            Ny = gridNy(iterNy)

            for iterP = 1 : length(gridP)

                p = gridP(iterP)

                %% prepare vectorized system
                Nx = Ny + 1;

                Nx0   = Nx * p;
                NyT   = Ny * T;
                NxT   = Nx * T;
                NxTp  = Nx * (T + p);

                kappa1 = .2^2;
                kappa2 = .5^2;
                kappa3 = 2;

                minnesotaPrior = ones(Ny,Ny);
                for i = 1 : Ny
                    for j = 1 : Ny
                        if i ~= j
                            minnesotaPrior(i,j) = kappa2;
                        end
                    end
                end
                minnesotaPrior = kappa1 .* minnesotaPrior;
                minnesotaPrior = minnesotaPrior ./ permute(1:p, [1 3 2]).^kappa3;
                agap           = minnesotaPrior .* randn(Ny,Ny,p);

                a                = zeros(Nx,Nx,p);
                a(1,1,1)         = 1;
                a(2:end,2:end,:) = agap;

                b = randn(Nx);
                b = chol(b * b')';

                c = [ones(Ny,1) eye(Ny, Ny)];

                % prior for x0; recall ordering is from x(-p+1) to x(0)
                X0                     = (1:Nx0)' * 100; % zeros(Nx0,1);
                cholsigX0              = 1e2 * eye(Nx0);
                cholsigX0(1:Nx:Nx0,1)  = 1e3;
                invcholsigX0           = inv(cholsigX0);

                % create 3D state space matrices
                aaa    = repmat(a, 1, 1, 1, T);
                ccc    = repmat(c, 1, 1, T);
                bbb    = repmat(b, [1 1 T]);
                invbbb = repmat(inv(b), [1 1 T]);

                %% construct matrix-form state space

                XX0   = sparse(1:Nx0, 1, X0, NxTp, 1);

                %% AA

                % build sequentially: first unit diagonal
                arows1     = 1 : NxTp;
                acols1     = 1 : NxTp;

                arows2     = repmat((1 : Nx)', 1, Nx * p);
                arows2     = Nx0 + arows2 + permute(Nx * (0 : T - 1), [1 3 2]);
                acols2     = repmat(1 : Nx * p, Nx,1) + permute(Nx * (0 : T - 1), [1 3 2]);

                arows      = [arows1(:); arows2(:)];
                acols      = [acols1(:); acols2(:)];

                values1    = ones(NxTp,1);
                values2    = -aaa(:,:,p:-1:1,:);
                values     = [values1(:); values2(:)];

                AA         = sparse(arows, acols, values);

                %% BB

                brows0  = repmat((1 : Nx0)', 1 , Nx0);
                brows1  = Nx0 + repmat((1 : Nx)', 1 , Nx) + permute(Nx * (0 : T-1), [1 3 2]);
                brows   = [brows0(:); brows1(:)];

                bcols0  = repmat((1 : Nx0), Nx0, 1);
                bcols1  = Nx0 + repmat((1 : Nx), Nx, 1) + permute(Nx * (0 : T-1), [1 3 2]);
                bcols   = [bcols0(:); bcols1(:)];

                bvalues = [cholsigX0(:); bbb(:)];

                BB      = sparse(brows, bcols, bvalues);

                invbvalues = [invcholsigX0(:); invbbb(:)];
                invBB      = sparse(brows, bcols, invbvalues);


                %% CC
                crows     = repmat((1 : Ny)', 1 , Nx, T) + permute(Ny * (0 : T-1), [1 3 2]);
                ccols     = Nx0 + repmat(1 : NxT, Ny, 1);
                CC        = sparse(crows(:), ccols(:), ccc(:), NyT, NxTp);

                %% simulate data
                % prepare
                rndStream  = getDefaultStream;
                xshocks    = randn(rndStream, NxTp, 1);
                X          = AA \ (XX0 + BB * xshocks(:));
                Y          = CC * X;
                
                YXprecsam0 = @() abcYXprecisionsampler(Y,XX0,AA,invBB,CC,rndStream);
                [~, QQ, RR1]  = YXprecsam0();
                YXprecsam = @() abcYXprecisionsampler(Y,XX0,AA,invBB,CC,rndStream,QQ,RR1);

                PS0yxtimes(iterP, iterNy, iterT) = timeit(YXprecsam0, 1);
                PSyxtimes(iterP, iterNy, iterT)  = timeit(YXprecsam, 1);


                y   = reshape(Y, Ny, T);
                aaa = aaa(:,:,p:-1:1);

                precsam0 = @() ALBCprecisionsampler(aaa,invbbb,ccc,y,X0,invcholsigX0,rndStream);
                [~, CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx]  = precsam0();
                precsam = @() ALBCprecisionsampler(aaa,invbbb,ccc,y,X0,invcholsigX0,rndStream,CC,QQ,RR1,arows, acols, asortndx, brows, bcols, bsortndx); 

                PS0times(iterP, iterNy, iterT) = timeit(precsam0, 1);
                PStimes(iterP, iterNy, iterT)  = timeit(precsam, 1);


                %% DK application
                Nw         = Nx;
                Nstates    = Nx * p;

                % construct index to remap prior
                ndx0 = (1:Nx0)';
                ndx0 = reshape(ndx0,Nx,p);
                ndx0 = fliplr(ndx0);
                ndx0 = ndx0(:);

                x0companion        = X0(ndx0);
                cholsigx0companion = cholsigX0(ndx0,ndx0);
                sigx0companion     = cholsigx0companion * cholsigx0companion';

                Acompanion         = zeros(Nstates);
                Acompanion(1:Nx,:) = reshape(a, Nx, Nstates);
                Acompanion(Nx+1:Nstates,1:Nx*(p-1)) = eye(Nx*(p-1));

                Bcompanion         = zeros(Nstates,Nw);
                Bcompanion(1:Nx,:) = b;

                Ccompanion         = zeros(Ny, Nstates);
                Ccompanion(:,1:Nx) = c;

                Acompanion = repmat(Acompanion, [1 1 T]);
                Bcompanion = repmat(Bcompanion, [1 1 T]);
                Ccompanion = repmat(Ccompanion, [1 1 T]);

                Ydata      = reshape(Y, Ny, T);
                dk         = @() abcDisturbanceSmoothingSampler1drawSLIM(Acompanion, Bcompanion, Ccompanion, Ydata, x0companion, cholsigx0companion, [], rndStream);

                DKtimes(iterP, iterNy, iterT) = timeit(dk, 1);

            end % p

        end % Ny

    end % T


    %% store results
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
    varlist = {'grid*', '*times', 'thisArch', 'thisVer', 'thisSys', 'thisBrand', ...
        'doSingleThread', 'usedThreads', 'availableThreads'};
    matname = sprintf('TRENDVARgapPStimes%sThreads%dof%d', thisBrand, usedThreads, availableThreads);
    save(matname, varlist{:});

    if quicky
        return
    end

end % doSingleThread
