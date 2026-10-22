%% ========================================================================
% TEST SUITE RUNNER
% ========================================================================

function nFail = runAllTests(varargin)
% Run the test tiers in order, log everything, and leave one report
%
% Inputs (name/value):
%   'Tier'  - highest tier to run, 0 to 4. Default 4.
%               0  environment: toolchain, data, shipped results   seconds
%               1  unit tests, no solver                           seconds
%               2  parity gate: realtimeSimulation vs the ref run   ~4 min
%               3  smoke tests: full year, and the tiled cascade   ~15 min
%               4  table reproduction, from the manifest            hours
%   'Only'  - tags for tier 4, overriding the manifest flags. Default {}.
%   'PlanOnly' - true runs tiers 0-3 and only PLANS tier 4.
%
% Outputs:
%   nFail  - total failed tiers; 0 means everything passed
%
% Usage:
%   runAllTests('Tier',1)                 % seconds, no solver needed
%   runAllTests('PlanOnly',true)          % everything cheap, plus the plan
%   runAllTests                           % the overnight run
%
% ORDER MATTERS. Tier 0 is a gate: if the environment is not ready nothing
% below it can mean anything, so a tier 0 failure stops the run. Tiers 1 to 4
% are independent, so a failure in one is recorded and the rest still run --
% an overnight pass should come back with the full picture, not with whatever
% happened to break first.
%
% Everything is written to logs/test_suite_<timestamp>.log as well as the
% console, because the point of an overnight run is reading it in the morning.

    p = inputParser;
    addParameter(p, 'Tier', 4);
    addParameter(p, 'Only', {});
    addParameter(p, 'PlanOnly', false);
    parse(p, varargin{:});
    opt = p.Results;

    ROOT = fileparts(fileparts(mfilename('fullpath')));
    addpath(ROOT); addpath(fullfile(ROOT,'functions')); addpath(fullfile(ROOT,'tests'));

    %% 1. Start the log
    logDir = fullfile(ROOT,'logs');
    if ~isfolder(logDir), mkdir(logDir); end
    logFile = fullfile(logDir, sprintf('test_suite_%s.log', ...
        datestr(now,'yyyymmdd_HHMMSS')));
    diary(logFile); diary on;

    t00 = tic;
    fprintf('\n########################################################\n');
    fprintf('# DDU CASCADE TEST SUITE\n');
    fprintf('# started  %s\n', datestr(now));
    fprintf('# root     %s\n', ROOT);
    if opt.PlanOnly
        fprintf('# tiers    0 to %d, tier 4 planned only\n', opt.Tier);
    else
        fprintf('# tiers    0 to %d\n', opt.Tier);
    end
    fprintf('# log      %s\n', logFile);
    fprintf('########################################################\n');

    names = {};
    fails = [];
    secs  = [];

    %% 2. Tier 0 is a gate
    t0 = tic;
    gateBroken = false;
    try
        n0 = checkEnvironment();
    catch err
        % The check itself threw, which is a bug in the check rather than a
        % verdict on the environment. Stopping here would cost an unattended
        % run the whole night over a faulty test, so record it and carry on.
        fprintf('\ncheckEnvironment ERRORED (a bug in the check, not a ');
        fprintf('verdict on this machine): %s\n', err.message);
        fprintf('Continuing anyway -- the later tiers will fail honestly if ');
        fprintf('the environment really is unusable.\n');
        n0 = 1;
        gateBroken = true;
    end
    names{end+1} = 'tier 0 environment';
    fails(end+1) = n0;
    secs(end+1)  = toc(t0);

    % Only a clean verdict of "not ready" stops the run.
    if n0 > 0 && ~gateBroken
        fprintf('\nSTOPPING: the environment is not ready, so no later tier ');
        fprintf('could be interpreted. Fix the items marked FAIL above.\n');
        nFail = finish(names, fails, secs, t00, logFile);
        return
    end

    %% 3. Tiers 1 to 4, each independent
    if opt.Tier >= 1
        [names, fails, secs] = tier(names, fails, secs, 'tier 1 unit tests', ...
            @() testUnits());
    end

    % verifyParity, the smoke tests and testSettlementNoise are SCRIPTS, so
    % they cannot be reached by a function handle. run() executes them, and
    % each goes through runScript so the `clear` at the top of a script wipes
    % only that helper's workspace and not the runner's bookkeeping.
    scripts = {};
    if opt.Tier >= 2
        scripts = [scripts; {'tier 2 parity gate', 'verifyParity.m'}];
    end
    if opt.Tier >= 3
        scripts = [scripts; {'tier 3 smoke, full year',    'smokeTestFullYear.m'}];
        scripts = [scripts; {'tier 3 smoke, tiled cascade','smokeTestScale.m'}];
        scripts = [scripts; {'tier 3 settlement noise',    'testSettlementNoise.m'}];
    end
    for q = 1:size(scripts,1)
        f = fullfile(ROOT,'tests',scripts{q,2});
        [names, fails, secs] = tier(names, fails, secs, scripts{q,1}, ...
            @() runScript(f));
    end

    if opt.Tier >= 4
        [names, fails, secs] = tier(names, fails, secs, 'tier 4 table reproduction', ...
            @() reproduceTables('Only',opt.Only,'PlanOnly',opt.PlanOnly));
    end

    nFail = finish(names, fails, secs, t00, logFile);
end


%% ----------------------------------------------------------------- helpers

function [names, fails, secs] = tier(names, fails, secs, name, fn)
% Run one tier, record its cost, and never let it take the suite down.
    fprintf('\n');
    t0 = tic;
    try
        n = fn();
        if isempty(n) || ~isnumeric(n), n = 0; end
    catch err
        fprintf('\n%s errored: %s\n', name, err.message);
        n = 1;
    end
    names{end+1} = name;
    fails(end+1) = n;
    secs(end+1)  = toc(t0);
end


function n = runScript(scriptPath)
% Execute a test SCRIPT. These signal failure by erroring, which tier()
% catches, so reaching the end means the test passed.
    run(scriptPath);
    n = 0;
end


function nFail = finish(names, fails, secs, t00, logFile)
    fprintf('\n########################################################\n');
    fprintf('# SUITE SUMMARY\n');
    fprintf('########################################################\n');
    fprintf('  %-34s %8s %10s\n', 'tier', 'failed', 'minutes');
    for q = 1:numel(names)
        if fails(q) == 0, mark = '  ok  '; else, mark = ' FAIL '; end
        fprintf('[%s] %-32s %8g %10.1f\n', mark, names{q}, fails(q), secs(q)/60);
    end
    nFail = sum(fails > 0);
    fprintf('\n%d tiers run, %d with failures, %.1f minutes total.\n', ...
        numel(names), nFail, toc(t00)/60);
    if nFail == 0
        fprintf('\nALL TIERS PASS.\n');
    else
        fprintf('\n%d TIER(S) FAILED -- see the detail above.\n', nFail);
    end
    fprintf('finished %s\n', datestr(now));
    fprintf('log: %s\n', logFile);
    diary off;
end
