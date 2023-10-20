function [x]=uamp_sbl_mmv(A,R)
[m1,n1]=size(A);
[~,L]=size(R);
Iter=200;
A=A/sqrt(sum(sum(abs(A).^2))/n1);
  
%% UAMP-SBL for MMV
[U,S,V]=svd(A);
lamdap=diag(S*S');
%% 参数初始化
s=zeros(m1,L);
iku=0.01;
ga=ones(n1,1);
r1=U'*R;A1=U'*A;
beta=1;
tx=ones(1,L);x=zeros(n1,L);
for it=1:Iter

beta_all=0;ga_all=0;
for l=1:L
    tp(:,l)=tx(:,l)*lamdap;
    p(:,l)=A1*x(:,l)-(tp(:,l).*s(:,l));
    vh(:,l)=tp(:,l)./(1+beta*tp(:,l));
    hhat(:,l)=((beta*tp(:,l)).*r1(:,l)+p(:,l))./(1+beta*tp(:,l));
    beta_all=beta_all+(sum((abs(r1(:,l)-hhat(:,l))).^2)+sum(vh(:,l)));
end
beta=(m1*L)/(beta_all);
for l=1:L
    ts(:,l)=1./(tp(:,l)+1/beta);
    s(:,l)=ts(:,l).*(r1(:,l)-p(:,l));
    tq(:,l)=lamdap'*ts(:,l)/(n1);
    tq(:,l)=1/(tq(:,l));
    q(:,l)=x(:,l)+tq(:,l)*(A1'*s(:,l));
    x(:,l)=q(:,l)./(1+tq(:,l).*ga);
    tx(:,l)=mean(tq(:,l)./(1+tq(:,l)*ga));
    %tx(:,l)=mean(tx(:,l));
    ga_all=ga_all+(abs(x(:,l)).^2+tx(:,l));
end
ga=(2*iku+1)./(ga_all/(L));
iku=0.5*sqrt(log(mean(ga))-mean(log(ga)));


end