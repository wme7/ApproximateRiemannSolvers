function res = FD_WENO5_EE1d(a,w,dx,fsplitMth)
% *************************************************************************
% Coded by Manuel A. Diaz, 02.10.2012, NTU Taiwan.
% Last update on 2016.04.29, NHRI Taiwan.
% *************************************************************************
%
% Domain cells (I{i}) reference:
%
%                |           |   u(i)    |           |
%                |  u(i-1)   |___________|           |
%                |___________|           |   u(i+1)  |
%                |           |           |___________|
%             ...|-----0-----|-----0-----|-----0-----|...
%                |    i-1    |     i     |    i+1    |
%                |-         +|-         +|-         +|
%              i-3/2       i-1/2       i+1/2       i+3/2
%
% ENO stencils (S{r}) reference:
%
%
%                               |___________S2__________|
%                               |                       |
%                       |___________S1__________|       |
%                       |                       |       |
%               |___________S0__________|       |       |
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                                      -|
%                                     i+1/2
%
%
%               |___________S0__________|
%               |                       |
%               |       |___________S1__________|
%               |       |                       |
%               |       |       |___________S2__________|
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                               |+
%                             i-1/2
%
% WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Careful!: by using circshift over our domain, we are implicitly creating a
% favorable code that automatically includes periodical boundary conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produce flux splitting 
switch fsplitMth
    case 'LF',  [v,u] = LF(a,w);    % Lax-Friedrichs (LF) Flux Splitting
    case 'RUS', [v,u] = Rusanov(w); % Rusanov (Rus) Flux Splitting
    case 'SHLL',[v,u] = SHLL(w);    % Split HLL (SHLL) flux 
    otherwise, error('Splitting method not set.');
end

%% Left flux reconstruction $f_{i+1/2}^{-}$
vmm = circshift(v,[0 2]);
vm  = circshift(v,[0 1]);
%v  = circshift(v,[0 0]);
vp  = circshift(v,[0 -1]);
vpp = circshift(v,[0 -2]);

% Polynomials
p0n = (2*vmm - 7*vm + 11*v)/6;
p1n = ( -vm  + 5*v  + 2*vp)/6;
p2n = (2*v   + 5*vp - vpp )/6;

% Smooth Indicators (Beta factors)
B0n = 13/12*(vmm-2*vm+v  ).^2 + 1/4*(vmm-4*vm+3*v).^2; 
B1n = 13/12*(vm -2*v +vp ).^2 + 1/4*(vm-vp).^2;
B2n = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2;

% Constants
d0n = 1/10; d1n = 6/10; d2n = 3/10; epsilon = 1e-6;

% Alpha weights 
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alpha2n = d2n./(epsilon + B2n).^2;
alphasumn = alpha0n + alpha1n + alpha2n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
w2n = alpha2n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n.*p0n + w1n.*p1n + w2n.*p2n;

%% Right flux reconstruction $f_{i+1/2}^{+}$ 
umm = circshift(u,[0 2]);
um  = circshift(u,[0 1]);
%u  = circshift(u,[0 0]);
up  = circshift(u,[0 -1]);
upp = circshift(u,[0 -2]);

% Polynomials
p0p = ( -umm + 5*um + 2*u  )/6;
p1p = ( 2*um + 5*u  - up   )/6;
p2p = (11*u  - 7*up + 2*upp)/6;

% Smooth Indicators (Beta factors)
B0p = 13/12*(umm-2*um+u  ).^2 + 1/4*(umm-4*um+3*u).^2; 
B1p = 13/12*(um -2*u +up ).^2 + 1/4*(um-up).^2;
B2p = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u -4*up+upp).^2;

% Constants
d0p = 3/10; d1p = 6/10; d2p = 1/10; epsilon = 1e-6;

% Alpha weights 
alpha0p = d0p./(epsilon + B0p).^2;
alpha1p = d1p./(epsilon + B1p).^2;
alpha2p = d2p./(epsilon + B2p).^2;
alphasump = alpha0p + alpha1p + alpha2p;

% ENO stencils weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;
w2p = alpha2p./alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p.*p0p + w1p.*p1p + w2p.*p2p;

%% Compute finite difference residual term, df/dx.
res = (hp-circshift(hp,[0,1])+hn-circshift(hn,[0,1]))/dx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flux splitting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lax-Friedrichs
function [Fp,Fm] = LF(a,q)
    global gamma
    
    % primary properties
    rho=q(1,:); Fm=q(2,:)./rho; E=q(3,:)./rho; 
    p=(gamma-1)*rho.*(E-0.5*Fm.^2);
    
    % flux vector of conserved properties
    F=[rho.*Fm; rho.*Fm.^2+p; Fm.*(rho.*E+p)];
    
    % Lax-Friedrichs flux
    Fp=          0.5*(F + a*q); 
    Fm=circshift(0.5*(F - a*q),[0,-1]); 
end

% Rusanov (or local Lax-Friedrichs)
function [Fp,Fm] = Rusanov(q)
    global gamma
    
    % primary properties
    rho=q(1,:); u=q(2,:)./rho; E=q(3,:)./rho; 
    p=(gamma-1)*rho.*(E-0.5*u.^2); a=sqrt(gamma*p./rho); 
    
    % flux vector of conserved properties
    F=[rho.*u; rho.*u.^2+p; u.*(rho.*E+p)];
    
    % positive and negative fluxes
    I=ones(3,1); % I = [1;1;1;] column vector
    Fp=          0.5*(F + I*a.*q); 
    Fm=circshift(0.5*(F - I*a.*q),[0,-1]); 
end

% Splitted HLL flux form Ref.[2]:
function [Fp,Fm] = SHLL(q)
    global gamma
    
    % primary properties
    rho=q(1,:); u=q(2,:)./rho; E=q(3,:)./rho; 
    p=(gamma-1)*rho.*(E-0.5*u.^2);
    
    % flux vector of conserved properties
    F=[rho.*u; rho.*u.^2+p; u.*(rho.*E+p)];
    
    % Mach number
    a=sqrt(gamma*p./rho); M = u./a; 
    
    % Produce corrections to Mach number
    M(M> 1)= 1; 
    M(M<-1)=-1;
    M2 = M.^2;
    
    % constant column vector [1;1;1]
    I = ones(3,1);
    
    Fp=           0.5*((I*(M+1)).*F + I*(a.*(1-M2)).*q); 
    Fm=circshift(-0.5*((I*(M-1)).*F + I*(a.*(1-M2)).*q),[0,-1]); 
end