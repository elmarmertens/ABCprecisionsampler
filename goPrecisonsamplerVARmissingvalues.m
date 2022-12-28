%% Apply precision-based ABC sampler to VAR(p) models with missing observations (and simulated data)


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

for missValShare = [.01 .05 .1 .2 .3 .5]
    
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

        [PStimes, PS0times, DKtimes, DKMisstimes]  = deal(NaN(length(gridP), length(gridNy), length(gridT)));

        for iterT = 1 : length(gridT)
            T = gridT(iterT)

            for iterNy = 1 : length(gridNy)
                Ny = gridNy(iterNy)


                for iterP = 1 : length(gridP)

                    p = gridP(iterP)

                    Nx = Ny;
                    Nw = Ny;

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
                    a = minnesotaPrior .* randn(Ny,Ny,p);


                    b = randn(Nx);
                    b = chol(b* b')';

                    c = eye(Ny, Nx);

                    % deterministic initial state
                    sig00        = zeros(Nx);
                    cholsig00    = zeros(Nx);
                    invcholsig00 = zeros(Nx);
                    x0           = randn(Nx,p);
                    xbar         = 1 * ones(Nx,1);


                    %% create 3D state space matrices
                    aaa = repmat(a, 1, 1, 1, T);
                    ccc = repmat(c, 1, 1, T);
                    bbb = repmat(b, [1 1 T]);
                    b0  = cholsig00;

                    %% companion form matrices

                    acompanion                       = zeros(Ny*p);
                    acompanion(1:Ny,:)               = reshape(a, Ny, Ny * p);
                    acompanion(Ny+1:end, 1:Ny*(p-1)) = eye(Ny*(p-1));

                    % abs(eig(acompanion))

                    acompanion1             = blkdiag(1, acompanion);
                    acompanion1(1+(1:Nx),1) = xbar;

                    bcompanion1             = zeros(Nx * p + 1, Nw);
                    bcompanion1(1+(1:Nx),:) = b;

                    x0companion1 = cat(1, 1, x0(:));
                    cholsig00companion1 = zeros(1 + Ny * p);

                    % expand all to 3D matrices
                    acompanion1 = repmat(acompanion1, [1 1 T]);
                    bcompanion1 = repmat(bcompanion1, [1 1 T]);

                    ccompanion1 = cat(2, zeros(Ny,1,T), ccc, zeros(Ny,Ny * (p-1),T));

                    
                    %% construct matrix-form state space
                    NyT   = Ny * T;
                    NwT   = Nw * T;
                    NxT   = Nx * T;
                    XX0   = repmat(xbar, T, 1);

                    % adjust for initial conditions
                    for k = 1 : p
                        xndx                       = (k-1) * Ny + (1 : Ny);
                        theseInitialLags           = p - k + 1;
                        thisA                      = reshape(a(:,:,k:p), Ny, Ny * theseInitialLags);
                        thisX0                     = reshape(x0(:,1:p-k+1), Ny * theseInitialLags, 1);
                        XX0(xndx)                  = XX0(xndx) + thisA * thisX0;
                    end

                    %% AA

                    % build sequentially: first unit diagonal
                    arows     = 1 : NxT;
                    acols     = 1 : NxT;
                    values    = ones(NxT,1);

                    % add p lags (sequentially)
                    for k = 1 : p
                        theserows   = repmat((1 : Nx)', 1 , Nx, T - k);
                        theserows   = theserows + permute(Nx * (k : T-1), [1 3 2]);
                        arows       = [arows(:); theserows(:)];

                        thesecols   = repmat(1 : Nx * (T - k), Nx, 1);
                        acols       = [acols(:); thesecols(:)];
                        values      = [values(:); -reshape(aaa(:,:,k,1:T-k), Nx * Nx * (T - k), 1)];
                    end

                    AA        = sparse(arows, acols, values);

                    %% BB
                    brows  = repmat((1 : Nx)', 1 , Nw, T);
                    brows  = brows + permute(Nw * (0 : T-1), [1 3 2]);
                    brows  = brows(:);

                    bcols   = repmat(1 : NwT, Nx, 1);
                    bcols   = bcols(:);

                    BB      = sparse(brows, bcols, bbb(:));

                    invbbb  = repmat(inv(b), [1 1 T]);



                    %% CC
                    crows     = repmat((1 : Ny)', 1 , Nx, T);
                    crows     = crows + permute(Ny * (0 : T-1), [1 3 2]);
                    ccols     = repmat(1 : NxT, Ny, 1);
                    CC        = sparse(crows(:), ccols(:), ccc(:), NyT, NxT);


                    %% simulate data

                    % prepare
                    rndStream  = getDefaultStream;

                    Nynan     = ceil(missValShare * NyT);
                    % draw missing value indicator
                    ndxNaN       = randi(rndStream, NyT, Nynan, 1);
                    yNaN         = false(NyT, 1);
                    yNaN(ndxNaN) = true;
                    yNaNndx      = reshape(yNaN, Ny, T);


                    % simulate data

                    % draw shocks
                    xshocks    = randn(rndStream, T * Nx,1);

                    %% simulate stacked system

                    % simulate
                    X       = AA \ (XX0 + BB * xshocks);
                    Y       = CC * X;
                    Y(yNaN) = NaN;

                    %% timeit comparison

                    Ydata          = reshape(Y, Ny, T);
                    Ydata(yNaNndx) = 0;

                    precsam0 = @() VARTVPSVprecisionsamplerNaN0const(aaa,invbbb,ccc,Ydata,yNaNndx,x0,xbar,rndStream);
                    [~, CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx]  = precsam0();

                    precsam  = @() VARTVPSVprecisionsamplerNaN0const(aaa,invbbb,ccc,Ydata,yNaNndx,x0,xbar,rndStream,CC,QQ,RR1,arows, acols, asortndx, brows, bcols, bsortndx);
                    dk       = @() abcDisturbanceSmoothingSamplerNaN1draw(acompanion1, bcompanion1, ccompanion1, Ydata, yNaNndx, x0companion1, cholsig00companion1, [], [], rndStream);
                    dkmiss   = @() stateABnanDraw1(acompanion1(:,:,1), bcompanion1(:,:,1), 1+(1:Ny), Ydata, yNaNndx, x0companion1, ones(Ny,T), rndStream);

                    PS0times(iterP, iterNy, iterT) =  timeit(precsam0, 1);
                    PStimes(iterP, iterNy, iterT)  =  timeit(precsam, 1);
                    DKtimes(iterP, iterNy, iterT)      =  timeit(dk, 1);
                    DKMisstimes(iterP, iterNy, iterT)  =  timeit(dkmiss, 1);

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
            'doSingleThread', 'usedThreads', 'availableThreads', 'missValShare'};
        matname = sprintf('VARmissingvaluesPStimes%sThreads%dof%dmissValShare%02d', thisBrand, usedThreads, availableThreads, floor(missValShare * 100));
        save(matname, varlist{:});

        if quicky
            return
        end

    end % doSingleThread
end % missValShare
