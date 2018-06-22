function dF = Fluctuations1d(qo,dx)

global gamma

% Compute boundary reconstructions: qR and qL
[ql,qr] = weno5(qo);
% [ql,qr] = muscl(qo,'MM');

% Compute fluctuations Amdq and Apdq
[~,~,amdq,apdq] = euler_hllc_1D(circshift(ql,[0,-1]),qr,gamma);
% [~,~,amdq,apdq] = euler_roe_1D(circshift(ql,[0,1]),qr,gamma);

% Compute total fluctuation
adq = totalFluctuation(ql,qr,gamma);

% Approximation of the spatial flux
dF = -(amdq + apdq + adq)/dx;
% dF = -(amdq + apdq)/dx;

end
