%% ========================================================================
% TEST TIER 0: ENVIRONMENT
% ========================================================================

function nFail = checkEnvironment()
% Everything a fresh clone needs, before any solver time is spent
%
% Outputs:
%   nFail  - number of failed checks; 0 means ready
%
% Usage:
%   checkEnvironment           % prints the table, errors if anything fails
%   n = checkEnvironment();    % returns the count instead of erroring
%
% Checks the toolchain, the input data, the shipped results and the path
% resolution, and reports EVERY problem in one pass rather than stopping at
% the first. This is the test that tells a new user why a clone will not run,
% instead of letting them find out 40 minutes into a sweep.

    ROOT = fileparts(fileparts(mfilename('fullpath')));
    addpath(ROOT); addpath(fullfile(ROOT,'functions'));
    addpath(fullfile(ROOT,'tests'));

    fprintf('\n===== TIER 0: ENVIRONMENT =====\n');
    fprintf('repository root: %s\n\n', ROOT);

    %% 1. Build the check list
    % Each row is {name, test, note, required}. required = false marks an
    % advisory check: something that limits what you can do but does not stop
    % the suite, so it is reported as WARN and does not fail the gate.
    spec = {};

    spec = [spec; {'MATLAB release', ...
        @() ~isempty(version('-release')), version('-release'), true}];
    spec = [spec; {'Parallel Computing Toolbox', ...
        @() ~isempty(ver('parallel')), ...
        'without it, pass Workers 1 to runSweep', false}];
    spec = [spec; {'gurobi on the path', ...
        @() exist('gurobi','file') > 0, 'add <gurobi>/matlab to the path', true}];
    spec = [spec; {'gurobi solves a trivial LP', ...
        @() gurobiWorks(), 'license may be missing or expired', true}];

    data = { ...
        'data/inputs.mat',              'water years and unit parameters'; ...
        'data/conditioned_shocks.mat',  'default realised shock path, one per year'; ...
        'data/archx_params.csv',        'fitted ARCH-X coefficients'; ...
        'data/bpa_pf24_lmp_wy2024.csv', 'the published BPA tariff series'};
    for q = 1:size(data,1)
        f = fullfile(ROOT, data{q,1});
        spec = [spec; {data{q,1}, @() isfile(f), data{q,2}, true}]; %#ok<AGROW>
    end

    % Every result the manifest says the paper depends on must ship the file
    % it is read from, plus the config needed to replay each of its tags. An
    % empty artefact means the result is reproduced on demand and nothing is
    % shipped for it, so there is nothing to check.
    M = paperTables();
    for q = 1:numel(M)
        if ~isempty(M(q).artefact)
            f = fullfile(ROOT, M(q).artefact);
            spec = [spec; {M(q).artefact, @() isfile(f), M(q).name, true}]; %#ok<AGROW>
        end
        for r = 1:numel(M(q).tags)
            c = fullfile(ROOT,'results','runs', ...
                         M(q).tags{r}, 'sweep_config.mat');
            spec = [spec; {sprintf('config for %s', M(q).tags{r}), ...
                @() isfile(c), 'reproduceTables replays this', true}]; %#ok<AGROW>
        end
    end

    bank = fullfile(ROOT,'data','scenario_bank_WY2024_N10_seed20260920_T8760.mat');
    spec = [spec; {'paired shock bank', @() isfile(bank), ...
        'every sweep draws its inflow paths from this', true}];

    ref = fullfile(ROOT,'tests','reference','s01_e050_DDU_k100_BON.mat');
    spec = [spec; {'parity reference job', @() isfile(ref), ...
        'verifyParity rebuilds this run and compares', true}];

    spec = [spec; {'python3 available', @() pythonWorks(), ...
        'only for figures/*.py; nothing else needs it', false}];
    spec = [spec; {'dataload rejects a bad root', @() dataloadGuards(), ...
        'a wrong path must raise dataload:noData', true}];

    %% 2. Run them
    n      = size(spec,1);
    passed = false(n,1);
    needed = cell2mat(spec(:,4));
    detail = cell(n,1);
    for q = 1:n
        try
            ok = spec{q,2}();
            passed(q) = ~isempty(ok) && all(logical(ok));
            detail{q} = spec{q,3};
        catch err
            passed(q) = false;
            detail{q} = err.message;
        end
        if passed(q)
            mark = '  ok  ';
        elseif needed(q)
            mark = ' FAIL ';
        else
            mark = ' warn ';
        end
        fprintf('[%s] %-42s %s\n', mark, spec{q,1}, detail{q});
    end

    %% 3. Verdict
    % Only a required check can fail the gate. An advisory one limits what you
    % can do -- render a Python figure, run a sweep in parallel -- without
    % making any result unobtainable, so it must not stop an overnight run.
    hardFail = ~passed &  needed;
    softFail = ~passed & ~needed;
    nFail    = sum(hardFail);
    fprintf('\n%d checks, %d failed, %d advisory.\n', n, nFail, sum(softFail));
    if any(softFail)
        fprintf('advisory: %s\n', strjoin(spec(softFail,1)', '; '));
    end
    if nFail > 0
        msg = sprintf('environment is not ready: %s', ...
                      strjoin(spec(hardFail,1)', '; '));
        if nargout == 0
            error('checkEnvironment:notReady', '%s', msg);
        end
        fprintf('%s\n', msg);
    else
        fprintf('TIER 0 PASS\n');
    end
end


function ok = gurobiWorks()
% One-variable LP: maximise x subject to x <= 1.
    m.A   = sparse(1,1,1,1,1);
    m.obj = 1;
    m.rhs = 1;
    m.sense = '<';
    m.modelsense = 'max';
    m.lb = 0;
    m.ub = 2;
    p.OutputFlag = 0;
    s = gurobi(m, p);
    ok = strcmp(s.status,'OPTIMAL') && abs(s.x - 1) < 1e-9;
end


function ok = pythonWorks()
    [st, ~] = system('python3 --version');
    ok = (st == 0);
end


function ok = dataloadGuards()
% dataload must assert on a root with no data/ folder. That assert is what
% turns a wrong path into a named error instead of a missing-field failure
% hundreds of lines later, so it is worth testing that it still fires.
    ok = false;
    try
        dataload(tempdir, 2024, 0.20, 7);
    catch err
        ok = strcmp(err.identifier, 'dataload:noData');
    end
end
