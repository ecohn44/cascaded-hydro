%% ========================================================================
% PARALLEL SWEEP DRIVER
% ========================================================================

function T = runSweep(varargin)
% Resumable parallel BON/SSH x DIU/DDU sweep, paper-consistent Sigma
%
% Inputs:
%   varargin  - name/value options; see the inputParser block below
%
% Outputs:
%   T         - the summary table written by collectResults
%
% Usage:
%   runSweep('Tag','N720_grid','MaxHours',720,'Scenarios',1:10, ...
%            'Epsilon',[0.01 0.05 0.10 0.20],'Kappa',[0.5 1 1.5],'Workers',7);
%
% Jobs: for each scenario and epsilon, DIU (kappa 0, which also serves as the
% DDU kappa=0 column; the control identity was verified in the smoke tests)
% plus DDU at each Kappa value, each solved with BON and SSH.
%
% Every finished job is saved to results/runs/<Tag>/jobs/ and is
% skipped on rerun, so the sweep can be stopped (Ctrl-C) and resumed. Full
% trajectories are kept only for SaveTrajEpsilon/SaveTrajKappa. At the end,
% collectResults(Tag) writes the summary CSV.

%% 1. Options

root=fileparts(mfilename('fullpath'));
addpath(root); addpath(fullfile(root,'functions'));
p=inputParser;
addParameter(p,'Tag','sweep');
addParameter(p,'MaxHours',720);
addParameter(p,'Scenarios',1:10);
addParameter(p,'Epsilon',[0.01 0.05 0.10 0.20]);
addParameter(p,'Kappa',[0.5 1 1.5]);
addParameter(p,'Methods',{'BON','SSH'});
addParameter(p,'Workers',7);
addParameter(p,'Horizon',24);
addParameter(p,'Theta',0);
addParameter(p,'TrackBand',0);
addParameter(p,'TrackSide','both');   % 'both' or 'below'
addParameter(p,'PriceMode','flat');   % 'flat', 'hourly', or 'bpa'
addParameter(p,'SigmaModel','binned'); % 'binned' (as-run) or 'archx'
addParameter(p,'PriceFlatten',0);     % 0 = TOU as-is, 1 = flat; compresses toward the mean
addParameter(p,'PriceNoiseSigma',0);  % log-sd of the real-time price departure (0 = off)
addParameter(p,'PriceNoiseRho',0.80); % hour-to-hour persistence of that departure
addParameter(p,'PriceNoiseSeed',[]);  % optional seed override (default: scenario seed)
addParameter(p,'ChanceMode','on');    % 'on' or 'off' (deterministic MPC benchmark)
addParameter(p,'WorldKappa',1);       % severity of the REALIZED world (1 = fitted DDU)
addParameter(p,'PolicyKappa',1);      % scale of the sigma the OPTIMIZER assumes (1 = fitted)
addParameter(p,'FloorPct',0.20);      % lower-envelope percentile (0.20 = P20)
addParameter(p,'PowerModel','fixed');   % 'fixed' or 'mccormick'
addParameter(p,'ChanceSides','lower');
addParameter(p,'Year',2024);
addParameter(p,'SaveTrajEpsilon',0.10);
addParameter(p,'SaveTrajKappa',[0 1]);
addParameter(p,'PlanOnly',false);   % true: return {cfg,pending} without running (used by run_batch)
addParameter(p,'BankFile',fullfile(root,'data','scenario_bank_WY2024_N10_seed20260920_T8760.mat'));
%% 2. Output folder and the sweep config
parse(p,varargin{:}); cfg=p.Results;

outDir=fullfile(root,'results','runs',cfg.Tag);
jobDir=fullfile(outDir,'jobs');
if ~isfolder(jobDir), mkdir(jobDir); end
save(fullfile(outDir,'sweep_config.mat'),'cfg');


% ---- job list
%% 3. Enumerate jobs, skipping the ones already on disk
jobs=struct('scen',{},'eps',{},'policy',{},'kappa',{},'method',{},'file',{});
for s=cfg.Scenarios(:)'
 for e=cfg.Epsilon(:)'
  pk={'DIU',0};
  for kk=cfg.Kappa(:)', pk(end+1,:)={'DDU',kk}; end %#ok<AGROW>
  for c=1:size(pk,1)
   for m=1:numel(cfg.Methods)
    f=sprintf('s%02d_e%03d_%s_k%03d_%s.mat',s,round(1000*e),pk{c,1}, ...
        round(100*pk{c,2}),cfg.Methods{m});
    jobs(end+1)=struct('scen',s,'eps',e,'policy',pk{c,1},'kappa',pk{c,2}, ...
        'method',cfg.Methods{m},'file',fullfile(jobDir,f)); %#ok<AGROW>
   end
  end
 end
end
done=arrayfun(@(j) isfile(j.file),jobs);
pending=jobs(~done);
% Run SSH jobs first (they are the long ones) so the pool stays balanced.
isSSH=strcmp({pending.method},'SSH');
pending=[pending(isSSH) pending(~isSSH)];
%% 4. Report the plan
fprintf('[%s] %s: %d jobs, %d done, %d pending, %d h each.\n', ...
    datestr(now,'HH:MM'),cfg.Tag,numel(jobs),nnz(done),numel(pending),cfg.MaxHours);

if cfg.PlanOnly
    T=struct('cfg',cfg,'pending',pending); return
end
C0=dataload(root,cfg.Year,cfg.FloorPct,7);
B=load(cfg.BankFile,'bank'); bank=B.bank;
if ~isempty(pending)
    if cfg.Workers>1
        pool=gcp('nocreate');
        if isempty(pool), parpool(cfg.Workers); end
        parfor q=1:numel(pending)
            sweepRunOne(pending(q),C0,bank,cfg);
        end
    else
        for q=1:numel(pending), sweepRunOne(pending(q),C0,bank,cfg); end
    end
end
%% 5. Collect
T=collectResults(cfg.Tag);
end
