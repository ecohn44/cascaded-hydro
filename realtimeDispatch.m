%% ========================================================================
% PER-HOUR DISPATCH LP
% ========================================================================

function [bonModel, B] = realtimeDispatch(t, S, ctx)
% Build the per-hour dispatch LP for the cascade
%
% Inputs:
%   t         - hour index within the water year, 1..T
%   S         - per-hour state: Vprev, uprev, Qprev, localPrev, lastResid2,
%               CrossingFirst, CrossingAny
%   ctx       - everything fixed for the whole run: the case C, the
%               configuration, the static constraint rows, the variable
%               index map and the sigma evaluator closure
%
% Outputs:
%   bonModel  - the Gurobi model under the Bonferroni margin
%   B         - every intermediate the solve phase needs: the constraint
%               triplets, the per-hour bounds, the head terms and the
%               diagnostics arrays. Returned whole rather than itemised, so
%               adding an intermediate here cannot silently fail to reach
%               the caller.
%
% Third of the three tiers, called ONCE PER HOUR. One rolling-horizon step:
% project the local inflow over the look-ahead, evaluate the forecast
% standard deviation under the configured uncertainty framework, shift the
% storage bounds by whatever margin the chance constraint requires, and
% assemble the Gurobi model.
%
% This function is PURE. It reads state and context and returns a model; it
% solves nothing, keeps no counters, and carries nothing between hours. The
% solve, the SSH handoff, the emergency recourse ladder and all state
% management live in realtimeSimulation.m, which is the only place the
% simulation advances.
%
% See realtime.md for the formulation and the variable ordering.

%% 1. Unpack
% Context and state are unpacked into local names identical to the ones the
% body uses, so the body below is byte-for-byte the code that ran as part of
% realtimeSimulation.m. Do not rename these.

Aeq = ctx.Aeq;
Aub0 = ctx.Aub0;
CrossingAny = S.CrossingAny;
Qhist = S.Qhist;
PrevPlanRel = S.PrevPlanRel;
uprev = S.uprev;
CrossingFirst = S.CrossingFirst;
H = ctx.H;
KAF = ctx.KAF;
Qprev = S.Qprev;
RHI = ctx.RHI;
RLO = ctx.RLO;
T = ctx.T;
Vprev = S.Vprev;
cfg = ctx.cfg;
lastResid2 = S.lastResid2;
lb = ctx.lb;
localPrev = S.localPrev;
lowerOnly = ctx.lowerOnly;
mcc = ctx.mcc;
n = ctx.n;
nnzRows = ctx.nnzRows;
nvar = ctx.nvar;
sigDIU = ctx.sigDIU;
sigEval = ctx.sigEval;
ub = ctx.ub;
z = ctx.z;


C = ctx.C;
CPOW = ctx.CPOW;
MC1 = ctx.MC1;
MC2 = ctx.MC2;
MC3 = ctx.MC3;
MC4 = ctx.MC4;
NR = ctx.NR;
POW = ctx.POW;
PROD = ctx.PROD;
RDN = ctx.RDN;
REL = ctx.REL;
RMIN = ctx.RMIN;
RNEG = ctx.RNEG;
RPOS = ctx.RPOS;
RPW = ctx.RPW;
RTN = ctx.RTN;
RTP = ctx.RTP;
RUP = ctx.RUP;
SHI = ctx.SHI;
SLO = ctx.SLO;
SPL = ctx.SPL;
TRK = ctx.TRK;
VOL = ctx.VOL;
ix = ctx.ix;
lmp = ctx.lmp;
maxRamp = ctx.maxRamp;
policy = ctx.policy;
trackBelow = ctx.trackBelow;

%% 2. Per-hour build
% Moved here verbatim from realtimeSimulation.m when the per-hour build was
% separated from the hour loop. The indentation is from that file and is kept
% deliberately: it is what made the move checkable line by line.
    t0=t-1;
    locHat=zeros(n,H);
    for i=1:n
        x=localPrev(i);
        for k=1:H
            x=C.ar(i,1)+C.ar(i,2)*x;
            locHat(i,k)=x;
        end
    end

    % Fix head from the implemented state before this rolling solve.
    headNow=zeros(1,n);
    for i=1:n
        p=C.sys(i);
        headNow(i)=max(p.fb0+p.fb_slope*Vprev(i) ...
            -(p.tw_a+p.tw_b*Qprev(i)),1);
    end

    lbNow=lb; ubNow=ub;
    headIntercept=zeros(1,n); mcEmergencyLower=zeros(n,H);
    if mcc
        for i=1:n
            p=C.sys(i);
            headIntercept(i)=p.fb0-(p.tw_a+p.tw_b*Qprev(i));
            if headIntercept(i)<1
                error('Head intercept below 1 ft at hour %d, unit %s.',t,C.units{i});
            end
        end
    end
    beq=zeros(n*H,1); bub=zeros(NR*n*H,1); obj=zeros(nvar,1);
    qLowH=zeros(H,n); qHighH=zeros(H,n);
    sdConst=zeros(n,H);
    pwI=zeros(6*n*H,1); pwJ=zeros(6*n*H,1); pwV=zeros(6*n*H,1); npw=0;
    rpI=zeros(2*maxRamp,1); rpJ=rpI; rpV=rpI; rpK=rpI; nrp=0;
    sdI=zeros(maxRamp,1); sdJ=sdI; sdV=sdI; nsd=0;
    for i=1:n
        p=C.sys(i);
        for k=1:H
            k0=k-1; tk=min(t+k0,T); r=(i-1)*H+k; b=(r-1)*NR;
            npw=npw+1; pwI(npw)=b+RPW; pwJ(npw)=ix(i,k,REL);
            if ~mcc
                pwV(npw)=-CPOW*headNow(i)/KAF;
            else
                % p <= CPOW*(headIntercept + fb_slope*V)*release/KAF,
                % with V*release replaced by the McCormick variable w.
                pwV(npw)=-CPOW*headIntercept(i)/KAF;
                npw=npw+1; pwI(npw)=b+RPW; pwJ(npw)=ix(i,k,PROD);
                pwV(npw)=-CPOW*p.fb_slope/KAF;
                rrDn=abs(p.rr_dn);
                rampLower=max(0,uprev(i)-k*rrDn)*KAF;
                rampUpper=min(p.umax,uprev(i)+k*p.rr_up)*KAF;
                yL=min(rampUpper,max(p.umin*KAF,rampLower)); yU=rampUpper;
                mcEmergencyLower(i,k)=rampLower;
                lbNow(ix(i,k,REL))=rampLower; ubNow(ix(i,k,REL))=rampUpper;
                ubNow(ix(i,k,PROD))=p.R*yU;
                npw=npw+1; pwI(npw)=b+MC1; pwJ(npw)=ix(i,k,VOL); pwV(npw)= yL;
                npw=npw+1; pwI(npw)=b+MC2; pwJ(npw)=ix(i,k,VOL); pwV(npw)= yU;
                npw=npw+1; pwI(npw)=b+MC3; pwJ(npw)=ix(i,k,VOL); pwV(npw)=-yL;
                npw=npw+1; pwI(npw)=b+MC4; pwJ(npw)=ix(i,k,VOL); pwV(npw)=-yU;
                bub(b+MC1)=0; bub(b+MC2)=p.R*yU; bub(b+MC3)=-p.R*yL; bub(b+MC4)=0;
            end

            j=C.upstream(i);
            if i==1
                routedKnown=C.routedExo(tk);
            elseif k0-C.lag(i)<0
                if C.lag(i)>0
                    routedKnown=known_out(j,t0+k0-C.lag(i));
                else
                    routedKnown=0;
                end
            else
                routedKnown=0;
            end
            beq(r)=(k==1)*Vprev(i)+(locHat(i,k)+routedKnown)*KAF;
            obj(ix(i,k,POW))=-lmp(tk);
            obj(ix(i,k,SPL))=C.spillPenalty;
            obj(ix(i,k,TRK))=cfg.Theta/p.R;
            obj(ix(i,k,SLO))=C.slackPenalty;
            obj(ix(i,k,SHI))=C.slackPenalty;

            rrDown=abs(p.rr_dn);
            bub(b+RUP)=p.rr_up*KAF+(k==1)*uprev(i)*KAF;
            bub(b+RDN)=rrDown*KAF-(k==1)*uprev(i)*KAF;
            bub(b+RPOS)=(k==1)*uprev(i)*KAF;
            bub(b+RNEG)=-(k==1)*uprev(i)*KAF;
            bub(b+RTP)=C.ref(tk,i)+cfg.TrackBand*p.R;
            if trackBelow, bub(b+RTP)=1e6; end
            bub(b+RTN)=-C.ref(tk,i)+cfg.TrackBand*p.R;
            bub(b+RMIN)=-p.umin*KAF;

            qLowH(k,i)=C.qLow(tk,i); qHighH(k,i)=C.qHigh(tk,i);
            scale=sqrt(k)/k; marginKnown=0;
            % Use the identical arithmetic path when kappa=0. This makes
            % the DDU=DIU control exact rather than relying on two
            % algebraically equivalent expressions in a degenerate LP.
            if strcmpi(policy,'DIU') || cfg.Kappa==0
                marginKnown=z*KAF*sqrt(k)*sigDIU(i);
                sdConst(i,k)=KAF*sqrt(k)*sigDIU(i);
            else
                s2prev=lastResid2(i);           % ARCH seed: last observed e^2
                for jj=0:k0
                    src=jj-C.lag(i);
                    if j==0 || src<0
                        if i==1
                            tt=min(t0+jj-C.lag(i),T-1);
                            du=abs(C.upReleaseExo(max(tt,0)+1) ...
                                -C.upReleaseExo(max(tt-1,0)+1));
                            uLev=C.upReleaseExo(max(tt,0)+1);
                        else
                            tt=t0+src;
                            du=abs(known_out(j,tt)-known_out(j,tt-1));
                            uLev=known_out(j,tt);
                        end
                        [sg,s2prev]=sigEval(i,uLev,du,s2prev);
                        marginKnown=marginKnown+z*KAF*scale*sg;
                        sdConst(i,k)=sdConst(i,k)+KAF*scale*sg;
                    else
                        % Planned upstream ramp inside the look-ahead:
                        % predicted from the previous rolling plan shifted
                        % one hour (fixed data at time t). Look-ahead index
                        % m=src+1 is real hour t+src; the previous plan
                        % stored hour t-1+kk in column kk.
                        mLook=src+1;
                        duPred=abs(plan_release(j,mLook)-plan_release(j,mLook-1));
                        [sg,s2prev]=sigEval(i,plan_release(j,mLook),duPred,s2prev);
                        marginKnown=marginKnown+z*KAF*scale*sg;
                        sdConst(i,k)=sdConst(i,k)+KAF*scale*sg;
                    end
                end
            end

            nominalMargin=z*KAF*sqrt(k)*sigDIU(i);
            if lowerOnly
                crossing=C.qLow(tk,i)+nominalMargin > C.qHigh(tk,i);
            else
                crossing=C.qLow(tk,i)+nominalMargin > C.qHigh(tk,i)-nominalMargin;
            end
            if crossing
                if k<=cfg.SSHLookaheadHours, CrossingAny(t)=true; end
                if k==1, CrossingFirst(t)=true; end
                bub(b+RLO)=-C.qLow(tk,i);
                bub(b+RHI)=C.qHigh(tk,i);
            elseif lowerOnly
                bub(b+RLO)=-C.qLow(tk,i)-marginKnown;
                bub(b+RHI)=C.qHigh(tk,i);
            else
                bub(b+RLO)=-C.qLow(tk,i)-marginKnown;
                bub(b+RHI)=C.qHigh(tk,i)-marginKnown;
            end
        end
    end

    % sparse() sums duplicate (row,col) pairs, matching the previous
    % Aub(r,c)=Aub(r,c)+coef accumulation.
    AubPow=Aub0+sparse(pwI(1:npw),pwJ(1:npw),pwV(1:npw),nnzRows,nvar);
    Aub=AubPow+sparse(rpI(1:nrp),rpJ(1:nrp),rpV(1:nrp),nnzRows,nvar);
    sdCoeff=sparse(sdI(1:nsd),sdJ(1:nsd),sdV(1:nsd),n*H,nvar);
    bonModel=makeModel(Aeq,Aub,beq,bub,obj,lbNow,ubNow);

%% 3. Pack
% Everything assigned above is returned. Over-returning is free;
% under-returning is a silent stale-value bug in the caller.
B.AubPow = AubPow;
B.beq = beq;
B.bub = bub;
B.headIntercept = headIntercept;
B.headNow = headNow;
B.lbNow = lbNow;
B.mcEmergencyLower = mcEmergencyLower;
B.nrp = nrp;
B.obj = obj;
B.qHighH = qHighH;
B.qLowH = qLowH;
B.rpI = rpI;
B.rpJ = rpJ;
B.rpK = rpK;
B.rpV = rpV;
B.sdCoeff = sdCoeff;
B.sdConst = sdConst;
B.ubNow = ubNow;

%% 4. Nested helpers
% Copies of the two nested helpers from realtimeSimulation.m. They are
% nested here, not passed as handles, so the call sites in section 2 stay
% byte-identical to the code that ran inside the pre-split file. They close over
% Qprev, PrevPlanRel, Qhist, H and T, which are unpacked in section 1.
    function value=plan_release(j,m)
        % Predicted release of unit j at look-ahead index m (real hour
        % t+m-1). m=0 is the implemented outflow Qprev. The previous plan
        % (made at hour t-1) stored hour t-1+kk in column kk, so hour
        % t+m-1 is column m+1; beyond the plan the last value persists.
        if m<=0
            value=Qprev(j);
        else
            value=PrevPlanRel(j,min(m+1,H));
        end
    end
    function value=known_out(j,tt0)
        tt0=min(tt0,T-1);
        if tt0>=0 && isfinite(Qhist(tt0+1,j))
            value=Qhist(tt0+1,j);
        else
            value=Qprev(j);
        end
    end
end
