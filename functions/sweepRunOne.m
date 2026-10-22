%% ========================================================================
% SINGLE SWEEP JOB
% ========================================================================

function sweepRunOne(job,C0,bank,cfg)
% Run and save one sweep job
t0=tic;
try
    C=C0; C.zseq=bank.zseq(:,:,job.scen);
    cr=struct('Horizon',cfg.Horizon,'Epsilon',job.eps,'Kappa',job.kappa, ...
        'Theta',cfg.Theta,'MaxHours',cfg.MaxHours,'ScenarioID',job.scen, ...
        'ScenarioSeed',bank.seeds(job.scen),'SolutionMethod',job.method, ...
        'SSHMaxIterations',20,'SSHPhiTolerance',5e-4, ...
        'SSHBoundaryTolerance',1e-5,'SSHLookaheadHours',cfg.Horizon, ...
        'SSHInteriorPoint','bonferroni','ProgressHours',1e9, ...
        'Verbose',false,'ChanceSides',cfg.ChanceSides,'TrackBand',cfg.TrackBand,'TrackSide',cfg.TrackSide,'PriceMode',cfg.PriceMode,'SigmaModel',cfg.SigmaModel,'PriceFlatten',cfg.PriceFlatten,'PriceNoiseSigma',cfg.PriceNoiseSigma,'PriceNoiseRho',cfg.PriceNoiseRho,'PriceNoiseSeed',cfg.PriceNoiseSeed,'ChanceMode',cfg.ChanceMode,'WorldKappa',cfg.WorldKappa,'PolicyKappa',cfg.PolicyKappa,'PowerModel',cfg.PowerModel);
    out=realtimeSimulation(C,job.policy,cr);
    tr=out.trajectory;
    small=struct('V',tr.V,'power',sum(tr.power,2), ...
        'joint_reliability',tr.joint_reliability, ...
        'bon_step_feasible',tr.bon_step_feasible,'ssh_step_feasible',tr.ssh_step_feasible, ...
        'bon_first_feasible',tr.bon_first_feasible,'ssh_first_feasible',tr.ssh_first_feasible, ...
        'head',mean(tr.head,1),'V0',min(max(C.obs(1,:),0),C.R));
    summary=out.summary; %#ok<NASGU>
    if abs(job.eps-cfg.SaveTrajEpsilon)<1e-9 && any(abs(job.kappa-cfg.SaveTrajKappa)<1e-9)
        traj=tr; %#ok<NASGU>
        save_job(job.file,summary,small,job,traj);
    else
        save_job(job.file,summary,small,job,[]);
    end
    fprintf('[%s] DONE s%d eps %.2f %s k%.1f %s: %.0f s\n',datestr(now,'HH:MM'), ...
        job.scen,job.eps,job.policy,job.kappa,job.method,toc(t0));
catch ME
    fid=fopen(strrep(job.file,'.mat','_error.txt'),'w');
    fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off')); fclose(fid);
    fprintf('[%s] FAILED s%d eps %.2f %s k%.1f %s: %s\n',datestr(now,'HH:MM'), ...
        job.scen,job.eps,job.policy,job.kappa,job.method,ME.message);
end
end

function save_job(f,summary,small,job,traj) %#ok<INUSD>
if isempty(traj)
    save(f,'summary','small','job');
else
    save(f,'summary','small','job','traj');
end
end
