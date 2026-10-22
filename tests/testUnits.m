%% ========================================================================
% TEST TIER 1: UNIT TESTS, NO SOLVER
% ========================================================================

function nFail = testUnits()
% Contracts and invariants that hold without solving anything
%
% Outputs:
%   nFail  - number of failed checks; 0 means all invariants hold
%
% Usage:
%   testUnits           % prints the table, errors if anything fails
%   n = testUnits();    % returns the count instead of erroring
%
% Runs in seconds and needs no Gurobi licence. Every check here is something
% that, if it broke, would produce a plausible-looking but wrong number
% rather than a crash -- which is the class of bug a results repository most
% needs to guard against.
%
% Each check is run inside its own try/catch so one failure does not hide
% the rest.

    ROOT = fileparts(fileparts(mfilename('fullpath')));
    addpath(ROOT); addpath(fullfile(ROOT,'functions'));
    addpath(fullfile(ROOT,'tests'));

    fprintf('\n===== TIER 1: UNIT TESTS =====\n\n');

    names = {};
    tests = {};

    names{end+1} = 'dataload honours its documented contract';
    tests{end+1} = @() t_dataloadContract(ROOT);

    names{end+1} = 'index map is a bijection onto 1..nvar';
    tests{end+1} = @() t_indexMap();

    names{end+1} = 'Bonferroni z matches the closed form';
    tests{end+1} = @() t_bonferroniZ();

    names{end+1} = 'forecastError respects its floor and ceiling';
    tests{end+1} = @() t_forecastClip();

    names{end+1} = 'forecastError at kappa = 0 drops decision dependence';
    tests{end+1} = @() t_forecastKappaZero();

    names{end+1} = 'sliceCase slices rows and leaves units alone';
    tests{end+1} = @() t_sliceCase(ROOT);

    names{end+1} = 'tileCase replicates units positionally';
    tests{end+1} = @() t_tileCase(ROOT);

    names{end+1} = 'scenarioMargin is deterministic in its seed';
    tests{end+1} = @() t_scenarioMargin();

    names{end+1} = 'sshProbabilityCache detects the independent case';
    tests{end+1} = @() t_probCacheIdentity();

    names{end+1} = 'every script resolves the repository root correctly';
    tests{end+1} = @() t_rootDepth(ROOT);

    names{end+1} = 'no script addpaths a folder that does not exist';
    tests{end+1} = @() t_addpathTargets(ROOT);

    n      = numel(tests);
    passed = false(n,1);
    detail = repmat({''}, n, 1);
    for q = 1:n
        try
            tests{q}();
            passed(q) = true;
        catch err
            passed(q) = false;
            detail{q} = err.message;
        end
        if passed(q), mark = '  ok  '; else, mark = ' FAIL '; end
        fprintf('[%s] %-52s %s\n', mark, names{q}, detail{q});
    end

    nFail = sum(~passed);
    fprintf('\n%d checks, %d failed.\n', n, nFail);
    if nFail > 0
        msg = sprintf('unit tests failed: %s', strjoin(names(~passed), '; '));
        if nargout == 0
            error('testUnits:failed', '%s', msg);
        end
        fprintf('%s\n', msg);
    else
        fprintf('TIER 1 PASS\n');
    end
end


%% ------------------------------------------------------------------ checks

function t_dataloadContract(ROOT)
% The shape and field names every downstream file assumes.
    C = dataload(ROOT, 2024, 0.20, 7);
    n = numel(C.units);
    assert(n == 4, 'expected 4 units, got %d', n);
    assert(numel(C.sys) == n, 'C.sys has %d elements, expected %d', numel(C.sys), n);

    required = {'umin','umax','rr_up','rr_dn','pmax','R','fb0','fb_slope','tw_a','tw_b'};
    have = fieldnames(C.sys);
    missing = setdiff(required, have');
    assert(isempty(missing), 'C.sys is missing %s', strjoin(missing, ', '));

    % rr_dn is stored NEGATIVE. Several call sites take abs() of it, and a
    % sign flip here silently disables the curtailment path.
    assert(all([C.sys.rr_dn] < 0), 'rr_dn must be negative for every unit');
    assert(all([C.sys.umax] > [C.sys.umin]), 'umax must exceed umin');

    for f = {'ref','qLow','qHigh','obs','local','outObs'}
        A = C.(f{1});
        assert(size(A,1) == 8760 && size(A,2) == n, ...
            '%s is %dx%d, expected 8760x%d', f{1}, size(A,1), size(A,2), n);
        assert(~any(isnan(A(:))), '%s contains NaN', f{1});
    end

    assert(all(C.qLow(:) >= -1e-9), 'qLow must be non-negative after clipping');
    assert(all(C.qHigh(:) - C.qLow(:) > -1e-9), 'qHigh must not fall below qLow');
    assert(all(C.ref(:) - C.qLow(:) > -1e-6), 'the reference must sit above the floor');
    assert(isequal(size(C.ar), [n 3]), 'C.ar must be %dx3', n);
end


function t_indexMap()
% ix(i,k,v) must hit every column of the stacked decision vector exactly
% once. An off-by-one here would overlay two variables and quietly corrupt
% the LP rather than raise anything.
    for mcc = [false true]
        n = 4; H = 24;
        if mcc, NV = 11; else, NV = 10; end
        ix = @(i,k,v) ((i-1)*H + (k-1))*NV + v;
        nvar = n*H*NV;
        hit = zeros(nvar,1);
        for i = 1:n
            for k = 1:H
                for v = 1:NV
                    c = ix(i,k,v);
                    assert(c >= 1 && c <= nvar, 'ix out of range: %d', c);
                    hit(c) = hit(c) + 1;
                end
            end
        end
        assert(all(hit == 1), ...
            'mcc=%d: %d columns hit more than once, %d never hit', ...
            mcc, sum(hit > 1), sum(hit == 0));
    end
end


function t_bonferroniZ()
% One lower tail per project at n = 4 and eps = 0.05 is the paper's
% controller, and 2.241403 is the value the parity gate reports.
    n = 4; eps = 0.05;
    zLower = -sqrt(2)*erfcinv(2*(1-eps/n));
    assert(abs(zLower - 2.241403) < 1e-6, ...
        'one-sided z is %.6f, expected 2.241403', zLower);
    zBoth = -sqrt(2)*erfcinv(2*(1-eps/(2*n)));
    assert(zBoth > zLower, 'splitting over 2n tails must raise the quantile');
    % Sanity against the standard normal: Phi(z) must equal 1 - eps/n.
    phi = 0.5*erfc(-zLower/sqrt(2));
    assert(abs(phi - (1 - eps/n)) < 1e-12, 'z does not invert the normal CDF');
end


function t_forecastClip()
% The floor and ceiling are what keep the ARCH recursion from running away.
    AX = fakeAX();
    omgK = [1 1 1 1];
    sfloor = 2.0; sceil = 3.0;
    % Huge drivers must still land on the ceiling.
    [sg, s2] = forecastError(1, 1e6, 1e6, 1e6, AX, omgK, 1, 1, sfloor, sceil);
    assert(abs(sg - sceil) < 1e-12, 'sg = %.6f, expected the ceiling %.6f', sg, sceil);
    assert(abs(s2 - sceil^2) < 1e-10, 's2 must be the CLIPPED variance');
    % Zero drivers must still land on the floor.
    [sg, s2] = forecastError(1, -1e6, -1e6, 0, AX, omgK, 1, 1, sfloor, sceil);
    assert(abs(sg - sfloor) < 1e-12, 'sg = %.6f, expected the floor %.6f', sg, sfloor);
    assert(abs(s2 - sfloor^2) < 1e-10, 's2 must be the CLIPPED variance');
end


function t_forecastKappaZero()
% kappa = 0 is the decision-INdependent control. The release level and the
% ramp must then make no difference at all; if they do, the DIU arm is not a
% control and every DDU-vs-DIU comparison in the paper is confounded.
    AX = fakeAX();
    omgK = [4 4 4 4];
    lo = 1e-6; hi = 1e6;      % wide floor/ceiling so nothing is clipped
    [a, ~] = forecastError(1,   0,   0, 0.5, AX, omgK, 0, 1, lo, hi);
    [b, ~] = forecastError(1, 150,  80, 0.5, AX, omgK, 0, 1, lo, hi);
    assert(abs(a - b) < 1e-12, ...
        'kappa = 0 still responds to the decision: %.12f vs %.12f', a, b);
    % With kappa = 1 the same change MUST move sigma, or the test above is
    % passing for the wrong reason.
    [c, ~] = forecastError(1,   0,   0, 0.5, AX, omgK, 1, 1, lo, hi);
    [d, ~] = forecastError(1, 150,  80, 0.5, AX, omgK, 1, 1, lo, hi);
    assert(d > c + 1e-9, 'kappa = 1 should raise sigma with release and ramp');
end


function t_sliceCase(ROOT)
    C  = dataload(ROOT, 2024, 0.20, 7);
    h0 = 4032; nH = 168;
    Cs = sliceCase(C, h0, nH);
    for f = {'ref','qLow','qHigh','obs','local','outObs'}
        assert(size(Cs.(f{1}),1) == nH, '%s has %d rows, expected %d', ...
            f{1}, size(Cs.(f{1}),1), nH);
        assert(size(Cs.(f{1}),2) == numel(C.units), '%s lost a unit', f{1});
        assert(isequal(Cs.(f{1})(1,:), C.(f{1})(h0,:)), ...
            '%s does not start at hour %d', f{1}, h0);
    end
    assert(numel(Cs.sys) == numel(C.sys), 'slicing must not touch the units');
    assert(Cs.sliceStart == h0 && Cs.sliceHours == nH);
end


function t_tileCase(ROOT)
% Tiling is what the scalability study uses. Since C.sys became a struct
% array, replica b must be a positional copy -- and that is exactly the
% property the old name-keyed form got right only by accident of naming.
    C  = dataload(ROOT, 2024, 0.20, 7);
    n0 = numel(C.units);
    Ct = tileCase(C, 3, 'LocalMode','zeromean', 'Quiet',true);
    assert(numel(Ct.units) == 3*n0, 'expected %d units, got %d', 3*n0, numel(Ct.units));
    assert(numel(Ct.sys)   == 3*n0, 'sys has %d elements, expected %d', ...
        numel(Ct.sys), 3*n0);
    for b = 1:3
        for i = 1:n0
            assert(isequal(Ct.sys((b-1)*n0+i), C.sys(i)), ...
                'replica %d unit %d does not match the original', b, i);
        end
    end
    % The cascade must be one chain, not three disconnected ones.
    assert(isequal(Ct.upstream, [0, 1:(3*n0-1)]), 'tiled topology is not a chain');
end


function t_scenarioMargin()
% The multiplier is a property of the sampling design, not of the inflow
% scenario, so two calls with the same seed must agree exactly. If they do
% not, the arms of the benchmark table are no longer internally paired.
    n = 4; R = eye(n);
    z1 = scenarioMargin(n, R, 200, 0, 20260928);
    z2 = scenarioMargin(n, R, 200, 0, 20260928);
    assert(z1 == z2, 'same seed gave %.12f then %.12f', z1, z2);
    z3 = scenarioMargin(n, R, 200, 0, 20260929);
    assert(z3 ~= z1, 'a different seed should give a different multiplier');
    assert(z1 > 0 && isfinite(z1), 'the multiplier must be a positive scalar');
end


function t_probCacheIdentity()
% With R = I the joint rectangle probability is an exact product of
% univariate normal probabilities, with closed-form gradients and no
% sampling. The cache must take that branch, because the SSH cuts and the
% kappa = 0 control both rely on repeated evaluations being bit-identical.
    c = sshProbabilityCache(eye(4));
    assert(c.identity, 'R = I was not recognised as the independent case');
    assert(c.n == 4, 'cache reports n = %d, expected 4', c.n);
    Rc = eye(4); Rc(1,2) = 0.3; Rc(2,1) = 0.3;
    c2 = sshProbabilityCache(Rc);
    assert(~c2.identity, 'a correlated R must not take the independent branch');
end


function t_rootDepth(ROOT)
% Guards a bug class that has already bitten this repository six times.
%
% A file that computes its root as fileparts(mfilename('fullpath')) is
% asking for its OWN folder; each extra fileparts goes one level up. So the
% number of fileparts calls must equal one plus the file's depth below the
% repository root. Get it wrong and the path still looks plausible --
% results land in functions/results/, or data is read from the repository's
% parent -- and nothing errors until much later, or never, because a stray
% copy of the data happens to sit there.
%
% Checked statically, by reading the files, because the alternative is
% running every driver.
    folders = {'', 'functions', 'tests', 'figures', 'reproduce'};
    offenders = {};
    for q = 1:numel(folders)
        d = fullfile(ROOT, folders{q});
        if ~isfolder(d), continue; end
        L = dir(fullfile(d, '*.m'));
        depth = numel(strsplit(folders{q}, filesep));
        if isempty(folders{q}), depth = 0; end
        for r = 1:numel(L)
            txt = fileread(fullfile(d, L(r).name));
            % Each assignment of the form  x = fileparts(...fileparts(mfilename
            hits = regexp(txt, ...
                '(\w+)\s*=\s*((?:fileparts\s*\(\s*)+)mfilename\(''fullpath''\)', ...
                'tokens');
            for h = 1:numel(hits)
                nfp = numel(strfind(hits{h}{2}, 'fileparts'));
                if nfp - 1 ~= depth
                    % figures/ and tests/ legitimately take their own folder
                    % first and go up on the next line; only flag a direct
                    % root assignment, which is named root/ROOT/projectRoot.
                    nm = hits{h}{1};
                    if any(strcmpi(nm, {'root','ROOT','projectRoot'}))
                        offenders{end+1} = sprintf('%s/%s: %s uses %d fileparts at depth %d', ...
                            folders{q}, L(r).name, nm, nfp, depth);   %#ok<AGROW>
                    end
                end
            end
        end
    end
    assert(isempty(offenders), 'wrong root depth:\n  %s', ...
        strjoin(offenders, sprintf('\n  ')));
end


function t_addpathTargets(ROOT)
% A stale addpath is silent: MATLAB warns at most, so a file that adds a
% folder that no longer exists keeps working only while the needed functions
% happen to be on the path already. This repository shipped five addpath
% calls to a src/ folder that does not exist in its layout.
    folders = {'', 'functions', 'tests', 'figures', 'reproduce'};
    offenders = {};
    for q = 1:numel(folders)
        d = fullfile(ROOT, folders{q});
        if ~isfolder(d), continue; end
        L = dir(fullfile(d, '*.m'));
        for r = 1:numel(L)
            txt = fileread(fullfile(d, L(r).name));
            want = regexp(txt, ...
                'addpath\s*\(\s*fullfile\s*\(\s*\w+\s*,\s*''([^'']+)''', ...
                'tokens');
            for h = 1:numel(want)
                sub = want{h}{1};
                if ~isfolder(fullfile(ROOT, sub))
                    offenders{end+1} = sprintf('%s/%s adds missing folder %s', ...
                        folders{q}, L(r).name, sub);                 %#ok<AGROW>
                end
            end
        end
    end
    assert(isempty(offenders), 'stale addpath:\n  %s', ...
        strjoin(offenders, sprintf('\n  ')));
end


%% ----------------------------------------------------------------- helpers

function AX = fakeAX()
% A minimal ARCH-X parameter set, built here rather than read from the CSV so
% the test states its own expected arithmetic and does not move when the fit
% is refreshed.
    AX.umin  = zeros(1,4);
    AX.urng  = 100*ones(1,4);
    AX.rmin  = zeros(1,4);
    AX.rrng  =  50*ones(1,4);
    AX.alpha = 0.25*ones(1,4);
    AX.gu    = 1.00*ones(1,4);
    AX.gr    = 2.00*ones(1,4);
    AX.sbar  = 9.00*ones(1,4);
end
