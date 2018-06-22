function [u,a,enthalpy,pl,pr] = roe_averages(q_l,q_r,gamma)

% Calculate Roe averages
rhoSqrtl = sqrt(q_l(1,:));
rhoSqrtr = sqrt(q_r(1,:));
rohsq2 = rhoSqrtl + rhoSqrtr;
pl = (gamma-1)*(q_l(3,:) - 0.5*(q_l(2,:).^2)./ q_l(1,:));
pr = (gamma-1)*(q_r(3,:) - 0.5*(q_r(2,:).^2)./ q_r(1,:));
u = (q_l(2,:)./rhoSqrtl + q_r(2,:)./rhoSqrtr)./ rohsq2;
enthalpy = ( (q_l(3,:)+pl)./rhoSqrtl + (q_r(3,:)+pr)./rhoSqrtr )./rohsq2;
a = sqrt((gamma-1)*(enthalpy - 0.5*u.^2));
