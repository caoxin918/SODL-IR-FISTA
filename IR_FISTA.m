function [xrf] = IR_FISTA(A, y, varargin)
%IR_FISTA 此处显示有关此函数的摘要
%   此处显示详细说明
[m,n]=size(A);

defaultLambda = []; % 默认 lambda 值
parser = inputParser;
addOptional(parser, 'lambda', defaultLambda);
parse(parser, varargin{:});
lambda = parser.Results.lambda;

if isempty(lambda)
    LambdaMax=max(abs(A'*y));
    para.lamb=1e-3*LambdaMax;
else
    para.lamb = lambda;
end
kt=@(x) (A'*(A*x-y));kt1=@(x)(abs(kt(x))-para.lamb);
kkt=@(x) (x~=0).*(kt(x)+para.lamb*sign(x))+kt1(x).*(kt1(x)>0).*(x==0);

%%

x=zeros(n,1);
err=1;
iter=1;
para.itermax=200;
para.ITERMAX=2000;
para.OutIterMax=500;
para.errtol=1e-6;
para.KKTErrtol=1e-10;
para.t=1;
para.zero_add=min(m,3000);
para.chose=1;
kktcond=1;
% L=max_svdnum(A,2)^2;%
% x=FISTA(A,y,para.lamb,L,1e-10,20);
%%  低精度解
time=tic;
Ax=A*x;
time1=[];
fvalue1=[];
len_W=[];
Liter=[];
Nzero=[];
while(kktcond>para.KKTErrtol&&iter<para.OutIterMax)
    gk=A'*(Ax-y);
    reswhole=gk+para.lamb*sign(x);
    W0=find(x==0&abs(gk)>(0.9999*para.lamb));
    [~,W0sort]=sort(abs(gk(W0)));
    W1=W0(W0sort(end-min(length(W0),para.zero_add)+1:end));
    W2=find(x~=0);
    W=union(W1,W2);
    kktcond=norm(reswhole(W2))+norm(abs(gk(x==0&abs(gk)>(para.lamb)))-para.lamb);
    AW=A(:,W);
    para.L=1*max_svdnum(AW,2)^2;
    para.itermax=min(para.ITERMAX,para.itermax*1.5);
    para.errtol=max(min(para.errtol/10,kktcond/10),para.KKTErrtol^2);
    %para.errtol=min(para.errtol/10,1e-11);
       
    if para.chose==1
        u =subFISTA(AW,y,x(W),para);
    else
        HWW=full(AW'*AW);r=AW'*(-y);
        u =subFISTA(HWW,r,x(W),para);
    end
    time1(end+1)=toc(time);
    fvalue1(end+1)=0.5*norm(AW*u-y)^2+para.lamb*sum(abs(u));
    Ax=Ax+AW*(u-x(W));
    
    % 测试soft函数
    %u = soft(u, lambda/L);

    x(W)=u;
    iter=iter+1;
    len_W(end+1)=length(W);
    Liter(end+1)=para.L;
    Nzero(end+1)=length(W2);
end
time=toc(time);
xrf=x;
kkt_xrf=norm(kkt(xrf))
f=@(x) (1/2*norm(A*x-y)^2+para.lamb*sum(abs(x)));
end
