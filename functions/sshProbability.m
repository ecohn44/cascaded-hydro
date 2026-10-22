%% ========================================================================
% JOINT RELIABILITY AND GRADIENTS
% ========================================================================

function [phi,gV,gSd,strain] = sshProbability(x,ctx,k)
% Joint storage reliability and exact local sensitivities
%
% The kth look-ahead storage vector is modeled as
%   V_real = V_nominal + error, error ~ N(0, D_k*R*D_k), D_k = diag(sd_k).
% ctx.sd_const and ctx.sd_coeff reproduce the accumulated uncertainty used
% by the Bonferroni constraints, so BON and SSH use the same model.
%
% k may be a vector of look-ahead hours. Then
%   phi    is numel(k) x 1
%   gV     is n x numel(k)   (d phi / d V_i)
%   gSd    is n x numel(k)   (d phi / d sd_i)
%   strain is n x numel(k)   (distance to nearest bound in sd units)
% For scalar k the shapes match the previous implementation (n x 1).
%
% Gradients are analytic. With standardized limits a_i=(qLow_i-V_i)/sd_i
% and b_i=(qHigh_i-V_i)/sd_i,
%   d phi/d V_i  = [f(a_i) C_i(a_i) - f(b_i) C_i(b_i)] / sd_i
%   d phi/d sd_i = [a_i f(a_i) C_i(a_i) - b_i f(b_i) C_i(b_i)] / sd_i
% where f is the standard normal density and C_i(c) is the probability
% that the other projects stay inside their limits given Z_i=c. With
% independent errors C_i is the product of the other marginal
% probabilities, so everything is exact and costs microseconds.
%
% Only the first output is computed when nargout==1.

if isfield(ctx,'prob') && ~isempty(ctx.prob)
    cache=ctx.prob;
else
    cache=sshProbabilityCache(ctx.Rcorr);
end

n=ctx.n;
k=k(:)';
m=numel(k);

V=reshape(x(ctx.v_index(:,k)),n,m);
rows=bsxfun(@plus,(1:n)',(k-1)*n);
sd=ctx.sd_const(:,k)+reshape(full(ctx.sd_coeff(rows(:),:)*x),n,m);
sd=max(sd,ctx.sd_floor);
qLow=ctx.q_low(k,:)';
qHigh=ctx.q_high(k,:)';

a=(qLow-V)./sd;
b=(qHigh-V)./sd;
invalid=any(qLow>=qHigh,1);            % 1 x m, empty rectangle

if cache.identity
    marginal=max(normal_cdf(b)-normal_cdf(a),0);
    phi=prod(marginal,1)';
else
    phi=sov_rectangle(a',b',cache.L,cache.W);
end
phi(invalid)=0;
phi=min(max(phi,0),1);

if nargout<2
    return
end

if cache.identity
    conditionalA=zeros(n,m);
    for i=1:n
        others=[1:i-1 i+1:n];
        if isempty(others)
            conditionalA(i,:)=1;
        else
            conditionalA(i,:)=prod(marginal(others,:),1);
        end
    end
    conditionalB=conditionalA;
else
    conditionalA=ones(n,m); conditionalB=ones(n,m);
    if n>1
        Wc=cache.W(:,1:max(n-2,0));
        for i=1:n
            ci=cache.cond(i); idx=ci.idx;
            ca=a(i,:)'; cb=b(i,:)';
            ca(~isfinite(ca))=0; cb(~isfinite(cb))=0;
            lowA=a(idx,:)'-ca*ci.r; highA=b(idx,:)'-ca*ci.r;
            lowB=a(idx,:)'-cb*ci.r; highB=b(idx,:)'-cb*ci.r;
            P=sov_rectangle([lowA;lowB],[highA;highB],ci.L,Wc);
            conditionalA(i,:)=P(1:m)';
            conditionalB(i,:)=P(m+1:end)';
        end
    end
end

fa=normal_pdf(a); fb=normal_pdf(b);
afa=a.*fa; bfb=b.*fb;
afa(~isfinite(afa))=0; bfb(~isfinite(bfb))=0;
gV=(fa.*conditionalA-fb.*conditionalB)./sd;
gSd=(afa.*conditionalA-bfb.*conditionalB)./sd;
gV(:,invalid)=0; gSd(:,invalid)=0;

strain=min(V-qLow,qHigh-V)./sd;
end

function p=normal_cdf(z)
p=0.5*erfc(-z/sqrt(2));
end

function f=normal_pdf(z)
f=exp(-0.5*z.^2)/sqrt(2*pi);
f(~isfinite(z))=0;
end

function p=sov_rectangle(A,B,L,W)
%SOV_RECTANGLE P(A<=Y<=B), Y~N(0,L*L'), rows of A/B are separate problems.
% Genz separation of variables with a fixed lattice W (N x d-1).
[m,d]=size(A);
if d==0
    p=ones(m,1);
    return
end
if d==1
    p=max(normal_cdf(B(:,1)/L(1,1))-normal_cdf(A(:,1)/L(1,1)),0);
    return
end
N=size(W,1);
Y=zeros(m,N,d-1);
f=ones(m,N);
s=zeros(m,N);
for kk=1:d
    if kk>1
        s=zeros(m,N);
        for j=1:kk-1
            s=s+L(kk,j)*Y(:,:,j);
        end
    end
    lo=normal_cdf(bsxfun(@minus,A(:,kk),s)/L(kk,kk));
    hi=normal_cdf(bsxfun(@minus,B(:,kk),s)/L(kk,kk));
    width=max(hi-lo,0);
    f=f.*width;
    if kk<d
        u=lo+bsxfun(@times,W(:,kk)',width);
        u=min(max(u,1e-300),1-1e-16);
        Y(:,:,kk)=-sqrt(2)*erfcinv(2*u);
    end
end
p=mean(f,2);
end
