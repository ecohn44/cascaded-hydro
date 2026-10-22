%% ========================================================================
% SSH INTERIOR POINT
% ========================================================================

function [slater,mode] = findSlater(slaterModel,gpar,ctx)
% Conservative feasible point for the SSH line search
%
% The SSH line search needs a point that is safely inside the joint chance
% set. A preselected minimum-release trajectory is sometimes too brittle:
% it can be infeasible even when the rolling LP has physically feasible
% decisions that sit well inside the storage chance region. We therefore
% first solve an auxiliary centering LP that maximizes the smallest
% Bonferroni-storage margin across all projects and look-ahead hours. If
% that fails, we fall back to the older minimum-release construction for
% diagnostics.

% Rows in realtimeSimulation.m are [balance; rolling inequalities]. Within each
% project-hour inequality block, rows 5 and 6 are storage lower/upper.
NR=10; if isfield(ctx,'NR'), NR=ctx.NR; end; RLO=5; RHI=6; Hf=ctx.Hfull;
storageRows=zeros(2*ctx.n*Hf,1);
q=0;
for i=1:ctx.n
    for k=1:Hf
        rr=(i-1)*Hf+k;
        q=q+1; storageRows(q)=ctx.n*Hf+(rr-1)*NR+RLO;
        q=q+1; storageRows(q)=ctx.n*Hf+(rr-1)*NR+RHI;
    end
end
if isfield(ctx,'lower_only') && ctx.lower_only
    storageRows=storageRows(1:2:end);   % margin to the P20 floor only
end

% Phase 1: find the most interior Bonferroni point. The added eta column
% tightens every storage chance row by eta kaf, while explicitly removing
% planned chance slack from these auxiliary rows. This keeps the Slater
% point conservative without pinning releases to an arbitrary path.
center=slaterModel;
etaIdx=numel(center.obj)+1;
center.A=[center.A sparse(size(center.A,1),1)];
center.obj=[zeros(numel(center.obj),1); -1];
center.lb=[center.lb; 0];
center.ub=[center.ub; inf];
if isfield(center,'varnames')
    center.varnames=[center.varnames; {'ssh_slater_margin_kaf'}];
end
% Build all centering rows in one block (row-by-row sparse growth is slow).
newRows=center.A(storageRows,1:etaIdx-1);
newRows(:,ctx.chance_slack_index)=0;
newRows=[newRows sparse(ones(numel(storageRows),1))];
nBefore=numel(center.rhs);
center.A=[center.A; newRows];
center.rhs=[center.rhs; center.rhs(storageRows)];
center.sense=[center.sense; repmat('<',numel(storageRows),1)];
if isfield(center,'constrnames')
    center.constrnames=[center.constrnames; arrayfun( ...
        @(q)sprintf('ssh_slater_center_%03d',q), ...
        (nBefore+1:nBefore+numel(storageRows))','UniformOutput',false)];
end

slater=gurobi(center,gpar);
centerStatus=slater.status;
if strcmp(slater.status,'OPTIMAL')
    eta=slater.x(etaIdx);
    slater.x=slater.x(1:ctx.nvar);
    mode=sprintf('max-centered Bonferroni point, eta %.4g kaf',eta);
    return
end

strict=slaterModel;
strict.lb(ctx.chance_slack_index)=0;
strict.ub(ctx.chance_slack_index)=0;
strict.lb(ctx.release_index)=ctx.slater_release;
strict.ub(ctx.release_index)=ctx.slater_release;
strict.lb(ctx.ramp_index)=ctx.slater_ramp;
strict.ub(ctx.ramp_index)=ctx.slater_ramp;

% First try the conservative release path with no spill.
candidate=strict;
candidate.lb(ctx.spill_index)=0;
candidate.ub(ctx.spill_index)=0;
slater=gurobi(candidate,gpar);
mode='minimum-release, zero-spill';
if strcmp(slater.status,'OPTIMAL'), return; end

% Spill is now permitted, but only because the minimum-release, zero-spill
% trajectory was infeasible. Minimize it directly to avoid gratuitous spill.
candidate=strict;
candidate.obj=zeros(size(candidate.obj));
candidate.obj(ctx.spill_index)=1;
slater=gurobi(candidate,gpar);
mode='minimum-release, required-spill';
if strcmp(slater.status,'OPTIMAL'), return; end

error(['No strict SSH Slater point exists under the conservative ' ...
    'centering or minimum-release constructions. Centering status: %s; ' ...
    'minimum-release status: %s.'],centerStatus,slater.status);
end
