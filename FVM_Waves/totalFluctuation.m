function adq = totalFluctuation(ql,qr,gamma)
%
%   "Internal" Riemann solver for the euler equations in 1D.
%   The riemann problem is solved by assuming a discontinuity at the
%   center of the i'th cell.
%
%   On input, ql contains the state vector at the left edge of each cell
%             qr contains the state vector at the right edge of each cell
%             auxl contains the auxiliary vector at the left edge of each cell
%             auxr contains the state vector at the right edge of each cell
%             maxnx is the number of physical points (without ghost cells)
%             mbc is the number of ghost cells
%             meqn is the number of equations
%             ixyz is the dimension index
%             mx is the size of the patch for the dimension corresponding
%               to the value of ixyz
%
%   On output, adq contains the decomposition of the flux difference
%              f(qr(i)) - f(ql(i)).
%
%   For the Euler equations q = (q1, q2, q3)^T, 
%   where q1 = rho, q2 = rho*u, q3 = E and
%              _                                              _  
%             |                       q2                       |
%   f(q) =    |        0.5(3-gamma)q2^2/q1 + (gamma-1)q3       |
%             |_  ( gamma*q3 - 0.5(gamma-1)q2^2/q1 ) * q2/q1  _|
%

mx = size(ql,2); adq=zeros(size(ql));

for i = 1:mx
    adq(1,i) = qr(2,i) - ql(2,i);
    adq(2,i) = 0.5*(3-gamma)*(qr(2,i)^2/qr(1,i) - ql(2,i)^2/ql(1,i)) + ...
               (gamma-1)*(qr(3,i) - ql(3,i));
    adq(3,i) = (gamma*qr(3,i)-0.5*(gamma-1)*qr(2,i)^2/qr(1,i))*qr(2,i)/qr(1,i) - ... 
               (gamma*ql(3,i)-0.5*(gamma-1)*ql(2,i)^2/ql(1,i))*ql(2,i)/ql(1,i);
end