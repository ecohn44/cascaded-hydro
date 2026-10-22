%% ========================================================================
% TEST TIER 4: TABLE REPRODUCTION
% ========================================================================

function nFail = reproduceTables(varargin)
% Rerun the paper's runs from scratch and diff them against what shipped
%
% Inputs (name/value):
%   'Only'     - cell array of tags to reproduce, overriding the manifest's
%                reproduce flags. Default {} uses the manifest.
%   'Tol'      - relative tolerance on the physical metrics. Default 1e-6.
%   'PlanOnly' - true prints the plan and the runtime estimate and runs
%                nothing. Use this before committing a night of solver time.
%   'Workers'  - parallel workers, overriding whatever the stored config
%                recorded on the machine that produced it. Default: as stored.
%
% Outputs:
%   nFail  - number of tags whose reproduction did not match
%
% Usage:
%   reproduceTables('PlanOnly',true)         % what would run, and for how long
%   reproduceTables                          % the manifest's default set
%   reproduceTables('Only',{'BPA_HL2_TH0000400'})
%
% HOW IT WORKS, AND WHY IT IS FAITHFUL. Each paper tag ships its own
% sweep_config.mat, which is the struct runSweep's own inputParser produced
% when that tag was run. Replaying those fields back through the same parser
% reproduces the run exactly, with no hand-transcribed argument list to drift
% out of date. Two fields are deliberately overridden:
%
%   Tag       -> REPRO_<tag>, so the shipped results are never overwritten.
%   BankFile  -> the repo-relative path, because the stored value is an
%                absolute path from the machine that produced the tag and
%                will not resolve anywhere else.
%
% WHAT COUNTS AS A MATCH. Same seeds and same solver should give the same
% numbers, but Gurobi is free to differ in the last bits, so physical metrics
% are compared at a relative tolerance and timings are excluded outright.
% Row identity (scenario, epsilon, policy, kappa, method) must match exactly;
% a mismatch there means the grid itself changed, which is a different and
% worse failure than a numerical one.
%
% RESUMABLE. runSweep skips finished jobs under the same Tag, so an
% interrupted pass continues where it stopped. Rerunning after a completed
% pass re-does only the comparison.

    p = inputParser;
    addParameter(p, 'Only', {});
    addParameter(p, 'Tol', 1e-6);
    addParameter(p, 'PlanOnly', false);
    addParameter(p, 'Workers', []);
    parse(p, varargin{:});
    opt = p.Results;

    ROOT = fileparts(fileparts(mfilename('fullpath')));
    addpath(ROOT); addpath(fullfile(ROOT,'functions'));
    addpath(fullfile(ROOT,'tests'));

    BANK = fullfile(ROOT,'data','scenario_bank_WY2024_N10_seed20260920_T8760.mat');
    % collectResults writes SHORTER timing names into the CSV than the job
    % summary struct carries, so a pattern written for verifyParity's field
    % names misses them. These are the only two timing columns in the CSV.
    SKIP_PAT = '^solve_mean_s$|^solve_p95_s$';
    KEY      = {'scen','eps','policy','kappa','method'};

    fprintf('\n===== TIER 4: TABLE REPRODUCTION =====\n\n');

    %% 1. Decide what to run
    M = paperTables();
    todo = {};
    for q = 1:numel(M)
        for r = 1:numel(M(q).tags)
            tag = M(q).tags{r};
            if strcmp(tag,'scalability'), continue; end
            if isempty(opt.Only)
                want = M(q).reproduce;
            else
                want = ismember(tag, opt.Only);
            end
            if want, todo{end+1} = tag; end     %#ok<AGROW>
        end
    end
    todo = unique(todo, 'stable');

    assert(~isempty(todo), 'reproduceTables:nothingToDo', ...
        'no tags selected. Check the manifest reproduce flags or the Only list.');
    assert(isfile(BANK), 'reproduceTables:noBank', ...
        'missing the shock bank: %s', BANK);

    hours = 0;
    fprintf('  %-28s %6s %8s\n', 'tag', 'jobs', 'hours');
    for q = 1:numel(todo)
        [jobs, hrs] = tagCost(M, todo{q});
        hours = hours + hrs;
        fprintf('  %-28s %6d %8.1f\n', todo{q}, jobs, hrs);
    end
    fprintf('  %-28s %6s %8.1f\n\n', 'TOTAL', '', hours);

    %% 2. Pre-flight: replay every tag's arguments without solving
    % runSweep's PlanOnly returns before it loads anything, so this exercises
    % the whole argument path -- the stored config, the field-to-name/value
    % rebuild, and runSweep's own inputParser -- in seconds. The point is that
    % an argument runSweep will not accept is found now, not six hours into an
    % unattended run that then produces nothing.
    fprintf('pre-flight: replaying each config through runSweep, no solving\n');
    dryFail = {};
    for q = 1:numel(todo)
        try
            runOneTag(ROOT, todo{q}, BANK, opt.Workers, true);
            fprintf('  [  ok  ] %s\n', todo{q});
        catch err
            dryFail{end+1} = sprintf('%s (%s)', todo{q}, err.message); %#ok<AGROW>
            fprintf('  [ FAIL ] %s: %s\n', todo{q}, err.message);
        end
    end
    if ~isempty(dryFail)
        msg = sprintf(['pre-flight failed, so nothing was solved: %s\n' ...
            'Fix this before committing solver time.'], strjoin(dryFail, '; '));
        if nargout == 0
            error('reproduceTables:preflight', '%s', msg);
        end
        fprintf('%s\n', msg);
        nFail = numel(dryFail);
        return
    end
    fprintf('pre-flight clean.\n');

    if opt.PlanOnly
        fprintf('\nPlanOnly: nothing was solved.\n');
        nFail = 0;
        return
    end

    %% 3. Rerun each tag into its own REPRO_ tag
    verdict = repmat({''}, numel(todo), 1);
    passed  = false(numel(todo),1);
    for q = 1:numel(todo)
        tag = todo{q};
        fprintf('\n----- %s -----\n', tag);
        try
            reproTag = runOneTag(ROOT, tag, BANK, opt.Workers, false);
            [passed(q), verdict{q}] = compareTag(ROOT, tag, reproTag, ...
                                                 opt.Tol, SKIP_PAT, KEY);
        catch err
            passed(q) = false;
            verdict{q} = err.message;
        end
        if passed(q), mark = '  ok  '; else, mark = ' FAIL '; end
        fprintf('[%s] %s: %s\n', mark, tag, verdict{q});
    end

    %% 4. Report
    fprintf('\n===== REPRODUCTION SUMMARY =====\n');
    for q = 1:numel(todo)
        if passed(q), mark = '  ok  '; else, mark = ' FAIL '; end
        fprintf('[%s] %-28s %s\n', mark, todo{q}, verdict{q});
    end
    nFail = sum(~passed);
    fprintf('\n%d tags, %d failed.\n', numel(todo), nFail);

    writeReport(ROOT, todo, passed, verdict);

    if nFail > 0
        msg = sprintf('reproduction failed for: %s', ...
                      strjoin(todo(~passed), ', '));
        if nargout == 0
            error('reproduceTables:mismatch', '%s', msg);
        end
        fprintf('%s\n', msg);
    else
        fprintf('TIER 4 PASS -- every selected table reproduces from a clean run.\n');
    end
end


%% ----------------------------------------------------------------- helpers

function [jobs, hrs] = tagCost(M, tag)
% The manifest costs a whole row, which may hold several tags, so split it.
    jobs = 0; hrs = 0;
    for q = 1:numel(M)
        if ismember(tag, M(q).tags)
            share = numel(M(q).tags);
            jobs = round(M(q).jobs / share);
            hrs  = M(q).hours / share;
            return
        end
    end
end


function reproTag = runOneTag(ROOT, tag, BANK, workers, planOnly)
% Replay a tag's own stored config into a REPRO_ tag.
% planOnly true enumerates the jobs and returns without solving.
    if nargin < 5, planOnly = false; end
    cfgFile = fullfile(ROOT,'results','runs',tag,'sweep_config.mat');
    assert(isfile(cfgFile), 'reproduceTables:noConfig', ...
        ['%s has no sweep_config.mat, so its settings are unknown. Copy it ' ...
         'from the machine that produced the tag.'], tag);
    S = load(cfgFile,'cfg');
    cfg = S.cfg;

    reproTag = ['REPRO_' tag];

    % The grid axes come from the SHIPPED CSV, not from the config.
    %
    % sweep_config.mat records the arguments of ONE runSweep call. A driver
    % that loops -- reproduce/kappaSweep calls runSweep once per scenario --
    % overwrites it every iteration, so the stored config describes only the
    % last call and replaying it reproduces a fraction of the tag. The CSV is
    % the ground truth for which jobs the tag actually contains.
    shipped = fullfile(ROOT,'results','runs',tag, ...
                       sprintf('summary_%s.csv', tag));
    if isfile(shipped)
        A = readtable(shipped);
        cfg.Scenarios = reshape(unique(A.scen), 1, []);
        cfg.Epsilon   = reshape(unique(A.eps), 1, []);
        cfg.Methods   = reshape(unique(cellstr(string(A.method))), 1, []);
        % DIU is enumerated automatically as kappa 0, so Kappa is the DDU set.
        kDDU = unique(A.kappa(strcmpi(string(A.policy),'DDU')));
        if ~isempty(kDDU)
            cfg.Kappa = reshape(kDDU, 1, []);
        end
    end

    % Rebuild the argument list from the parsed struct. Every field came out
    % of runSweep's own parser, so every field is a name it accepts.
    cfg.Tag      = reproTag;
    cfg.BankFile = BANK;
    if ~isempty(workers), cfg.Workers = workers; end

    cfg.PlanOnly = planOnly;

    fn = fieldnames(cfg);
    args = cell(1, 2*numel(fn));
    for k = 1:numel(fn)
        args{2*k-1} = fn{k};
        args{2*k}   = cfg.(fn{k});
    end

    if ~planOnly
        fprintf('  replaying %d settings into %s\n', numel(fn), reproTag);
    end
    runSweep(args{:});
end


function [ok, verdict] = compareTag(ROOT, tag, reproTag, tol, skipPat, key)
% Diff the regenerated summary against the shipped one, row by row.
    f1 = fullfile(ROOT,'results','runs',tag, ...
                  sprintf('summary_%s.csv', tag));
    f2 = fullfile(ROOT,'results','runs',reproTag, ...
                  sprintf('summary_%s.csv', reproTag));
    assert(isfile(f1), 'shipped summary missing: %s', f1);
    assert(isfile(f2), 'regenerated summary missing: %s', f2);

    A = readtable(f1);
    B = readtable(f2);

    % Row identity first. If the grids differ, the numbers are not comparable
    % and reporting a tolerance failure would be misleading.
    assert(height(A) == height(B), ...
        'row count %d -> %d: the grid changed, not just the numbers', ...
        height(A), height(B));
    A = sortrows(A, key);
    B = sortrows(B, key);
    for k = 1:numel(key)
        c = key{k};
        if iscell(A.(c)) || isstring(A.(c))
            same = all(strcmpi(string(A.(c)), string(B.(c))));
        else
            same = all(abs(double(A.(c)) - double(B.(c))) < 1e-12);
        end
        assert(same, 'the %s column does not line up after sorting', c);
    end

    % Then the numbers.
    cols = intersect(A.Properties.VariableNames, B.Properties.VariableNames);
    worst = 0; worstCol = ''; bad = {}; nChecked = 0; skipped = {};
    for k = 1:numel(cols)
        c = cols{k};
        if ~isempty(regexpi(c, skipPat, 'once'))
            skipped{end+1} = c;                  %#ok<AGROW>
            continue
        end
        a = A.(c); b = B.(c);
        if ~isnumeric(a) || ~isnumeric(b), continue; end
        nChecked = nChecked + 1;
        finite = isfinite(a) & isfinite(b);
        assert(all(isnan(a(~finite)) == isnan(b(~finite))), ...
            'column %s: NaN pattern differs', c);
        if ~any(finite), continue; end
        rel = max(abs(a(finite) - b(finite)) ./ max(abs(a(finite)), 1));
        if rel > worst, worst = rel; worstCol = c; end
        if rel > tol
            bad{end+1} = sprintf('%s (%.2e)', c, rel);   %#ok<AGROW>
        end
    end

    ok = isempty(bad);
    if ok && worst == 0
        % Bit-identical, which is the right expectation: one shock bank at a
        % fixed seed, Threads = 1, a fixed Gurobi seed, and a pure LP at every
        % step. Merely "close" would mean something non-deterministic crept in.
        verdict = sprintf('%d rows, %d columns BIT-IDENTICAL', ...
            height(A), nChecked);
    elseif ok
        verdict = sprintf('%d rows, %d columns match; worst rel diff %.2e on %s', ...
            height(A), nChecked, worst, worstCol);
    else
        verdict = sprintf('%d of %d columns exceed %.1e: %s', ...
            numel(bad), nChecked, tol, strjoin(bad, ', '));
    end
    fprintf('  %d rows, %d columns compared, %d excluded as timing\n', ...
        height(A), nChecked, numel(skipped));
end


function writeReport(ROOT, tags, passed, verdict)
% One CSV under results/tests/, so an overnight run leaves a readable artefact
% even if nobody watched the transcript, without putting test output among the
% published summaries.
    outDir = fullfile(ROOT,'results','tests');
    if ~isfolder(outDir), mkdir(outDir); end
    T = table(string(tags(:)), passed(:), string(verdict(:)), ...
              'VariableNames', {'tag','passed','detail'});
    f = fullfile(outDir, 'reproduction_report.csv');
    writetable(T, f);
    fprintf('\nreport written: %s\n', f);
end
