%% ========================================================================
% JOINT RELIABILITY CACHE
% ========================================================================

function cache = sshProbabilityCache(Rcorr, numPoints)
% Precompute everything sshProbability needs once
%
% Build this once per realtimeSimulation call (not once per rolling hour). It
% holds the Cholesky factors and the fixed quasi-Monte Carlo lattice used
% by the general correlated evaluator, and flags the independent case.
%
% Independent case (Rcorr = I, which is what realtimeSimulation currently uses
% because dataload does not set C.Rcorr): the joint rectangle
% probability is an exact product of univariate normal probabilities and
% its gradients are closed form. No sampling is involved.
%
% Correlated case: Genz separation-of-variables with a FIXED Richtmyer
% lattice plus antithetic points. Because the lattice never changes, the
% estimate is a deterministic, smooth function of the limits, so the SSH
% bisection converges and repeated evaluations are bit-identical (which
% the kappa-zero control requires). mvncdf in d>=4 is randomized QMC and
% has neither property.

if nargin<2 || isempty(numPoints), numPoints=8192; end
n=size(Rcorr,1);
Rcorr=(Rcorr+Rcorr')/2;
offDiagonal=Rcorr-diag(diag(Rcorr));

cache.n=n;
cache.R=Rcorr;
cache.identity=max(abs(offDiagonal(:)))<1e-12 && ...
    max(abs(diag(Rcorr)-1))<1e-12;
cache.L=[];
cache.W=zeros(0,max(n-1,0));
cache.cond=struct('idx',{},'r',{},'L',{});
if cache.identity
    return
end

cache.L=chol_psd(Rcorr);

primesList=[2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71];
dims=max(n-1,1);
if dims>numel(primesList)
    error('sshProbabilityCache supports at most %d projects.',numel(primesList)+1);
end
half=ceil(numPoints/2);
j=(1:half)';
W=mod(j*sqrt(primesList(1:dims)),1);
cache.W=[W;1-W];

for i=1:n
    idx=[1:i-1 i+1:n];
    r=Rcorr(idx,i);
    S=Rcorr(idx,idx)-r*r';
    cache.cond(i).idx=idx;
    cache.cond(i).r=r(:)';
    if isempty(idx)
        cache.cond(i).L=zeros(0,0);
    else
        cache.cond(i).L=chol_psd(S);
    end
end
end

function L=chol_psd(S)
S=(S+S')/2;
[L,p]=chol(S,'lower');
jitter=1e-12*max(1,max(abs(diag(S))));
while p>0
    [L,p]=chol(S+jitter*eye(size(S)),'lower');
    jitter=10*jitter;
    if jitter>1e-2
        error('Correlation matrix is not positive semidefinite.');
    end
end
end
