function [IIR] = bmt_fun(fir, k)
%BMT_FUN Given an FIR and desired truncation point k, this function spits
% out the reduced state-space IIR

[V, ~, b] = bmt_eigs(fir);

n = length(b)-1;

% state-space reduction
Ak = V(2:n,1:k)'*V(1:n-1,1:k);
Bk = V(1,1:k)';
Ck = b(2:end)*V(1:n,1:k);
Dk = b(1);

% reduced system state-space IIR
[bk,ak] = ss2tf(Ak,Bk,Ck,Dk);
IIR = dfilt.df1(bk,ak);

end

