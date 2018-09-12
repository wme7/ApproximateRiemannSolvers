function res = FV_WENO_EE1d(q,nx,dx,fluxMethod,Recon)
% Compute RHS of the semi-discrete form of the Euler equations.
global gamma

%   Flux at j+1/2
% 
%     j+1/2    Cell's grid: (assuming WENO5, R=3)
%   |   |   |
%   | wL|   |               {x=0}                     {x=L}
%   |  /|wR |           1   2 | 3   4   5        N-3 N-2|N-1  N
%   | / |\  |         |-o-|-o-|-o-|-o-|-o-| ... |-o-|-o-|-o-|-o-|
%   |/  | \ |             1   2   3   4            N-3 N-2 N-1  
%   |   |  \|
%   |   |   |       NC: Here cells 1 to 2 and N-1 to N are ghost cells
%     j  j+1            faces 2 and N-2, are the real boundary faces.
%
%   q = cat(3, r, ru, rv, E);
%   F = cat(3, ru, ru^2+p, ruv, ru(E+p));
%   G = cat(3, rv, ruv, rv^2+p, rv(E+p));

% Identify number of gost cells
switch Recon
    case {'WENO5','Poly5'}, R=3; % R: stencil size and number of gost cells
    case {'WENO7','Poly7'}, R=4;
    otherwise, error('reconstruction not available ;P');
end

%% 0. Set boundary conditions for Reiomann Problems: outflux BCs
for i=1:R-1
    q(:,i,:)=q(:,R,:); q(:,nx+1-i,:)=q(:,nx+1-R,:);	% Neumann BCs
end

%% 1. Reconstruct interface values: qL=q_{i+1/2}^{-} and qR=q_{i-1/2}^{+}

% Compute primitive variables at solution points
w(1,:) = q(1,:);
w(2,:) = q(2,:)./q(1,:);
w(3,:) = (gamma-1)*( q(3,:) - 0.5*(q(2,:).^2)./q(1,:));

% Produce reconstruction
switch Recon
    case 'WENO5', [wL,wR] = WENO5recon(w,nx);
    case 'WENO7', [wL,wR] = WENO7recon(w,nx);
    case 'Poly5', [wL,wR] = POLY5recon(w,nx);
    case 'Poly7', [wL,wR] = POLY7recon(w,nx);
    otherwise, error('reconstruction not available ;P');
end

% Compute conservative variables at faces
qL(1,:) = wL(1,:);
qL(2,:) = wL(2,:).*wL(1,:);
qL(3,:) = wL(3,:)./(gamma-1) + 0.5*wL(1,:).*wL(2,:).^2;

qR(1,:) = wR(1,:);
qR(2,:) = wR(2,:).*wR(1,:);
qR(3,:) = wR(3,:)./(gamma-1) + 0.5*wR(1,:).*wR(2,:).^2;

%% Compute finite volume residual term, df/dx.
res=zeros(size(w)); flux=zeros(size(qR));
for j = R:(nx-R), J=j+1-R; % for all interior faces 
    % compute flux at i+1/2
    %flux(:,j) = LFflux(qn(:,j-2),qp(:,j-2),gamma,smax);
    switch fluxMethod
        case 'ROE', flux(:,J) = ROEflux(qL(:,J),qR(:,J)); % Roe
        case 'RUS', flux(:,J) = RUSflux(qL(:,J),qR(:,J));  % Rusanov
        case 'HLLE',flux(:,J) = HLLEflux(qL(:,J),qR(:,J)); % HLLE
        case 'AUSM',flux(:,J) = AUSMflux(qL(:,J),qR(:,J)); % AUSM
        case 'HLLC',flux(:,J) = HLLCflux(qL(:,J),qR(:,J)); % HLLC
    end
    % Flux contribution to the residual of every cell
    res(:, j ) = res(:, j ) + flux(:,J)/dx;
    res(:,j+1) = res(:,j+1) - flux(:,J)/dx;
end

% Flux contribution of the LEFT MOST FACE: left face of cell j=1.
%flux(:,3) = LFflux(qn(:,R),qp(:,R),gamma,smax);
switch fluxMethod
    case 'ROE', flux(:,1) = ROEflux(qL(:,1),qR(:,1)); % Roe
    case 'RUS', flux(:,1) = RUSflux(qL(:,1),qR(:,1)); % Rusanov
    case 'HLLE',flux(:,1) = HLLEflux(qL(:,1),qR(:,1)); % HLLE
    case 'AUSM',flux(:,1) = AUSMflux(qL(:,1),qR(:,1)); % AUSM
    case 'HLLC',flux(:,1) = HLLCflux(qL(:,1),qR(:,1)); % HLLC
end
res(:,R) = res(:,R) - flux(:,1)/dx;
 
% Flux contribution of the RIGHT MOST FACE: right face of cell j=nx-1.
%flux(:,nx-2) = LFflux(qn(:,nx-1-2*R),qp(:,nx-1-2*R),gamma,smax);
J=nx+1-2*R;
switch fluxMethod
    case 'ROE', flux(:,J) = ROEflux(qL(:,J),qR(:,J)); % Roe
    case 'RUS', flux(:,J) = RUSflux(qL(:,J),qR(:,J)); % Rusanov
    case 'HLLE',flux(:,J) = HLLEflux(qL(:,J),qR(:,J)); % HLLE
    case 'AUSM',flux(:,J) = AUSMflux(qL(:,J),qR(:,J)); % AUSM
    case 'HLLC',flux(:,J) = HLLCflux(qL(:,J),qR(:,J)); % HLLC
end
res(:,nx+1-R) = res(:,nx+1-R) + flux(:,J)/dx;

end % FVM WENO

function Roe = ROEflux(qL,qR)
    % Roe flux function
    global gamma
    
    % Left state
    rL = qL(1);
    uL = qL(2)./rL;
    EL = qL(3)./rL;
    pL = (gamma-1)*( qL(3) - rL*uL*uL/2 );
    aL = sqrt(gamma*pL/rL);
    HL = ( qL(3) + pL )./rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)./rR;
    ER = qR(3)./rR;
    pR = (gamma-1)*( qR(3) - rR*uR*uR/2 );
    aR = sqrt(gamma*pR/rR);
    HR = ( qR(3) + pR )./rR;
    
    % First compute the Roe Averages
    RT = sqrt(rR/rL);
    r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    H = (HL+RT*HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-u*u/2) );
    
    % Differences in primitive variables.
    dr = rR - rL;
    du = uR - uL;
    dP = pR - pL;
    
    % Wave strength (Characteristic Variables).
    dV = [(dP-r*a*du)/(2*a^2); -( dP/(a^2)-dr); (dP+r*a*du)/(2*a^2)];
    
    % Absolute values of the wave speeds (Eigenvalues)
    ws = [ abs(u-a); abs( u ); abs(u+a) ];

    % Harten's Entropy Fix JCP(1983), 49, pp357-393.
    % There are various ways to implement the entropy fix. This is just one
    % example. Try turn this off. The solution may be more accurate.
    Da = max(0,4*((uR-aR)-(uL-aL))); if (ws(1)<Da/2); ws(1)=ws(1)*ws(1)/Da+Da/4; end
    Da = max(0,4*((uR+aR)-(uL+aL))); if (ws(3)<Da/2); ws(3)=ws(3)*ws(3)/Da+Da/4; end

    % Right eigenvectors
    R = [  1  ,  1  ,  1  ;
         u-a ,  u  , u+a ;
        H-u*a,u^2/2,H+u*a];
   
    % Compute the average flux.
    FL=[rL.*uL; rL.*uL.^2+pL; uL.*(rL.*EL+pL)];
    FR=[rR.*uR; rR.*uR.^2+pR; uR.*(rR.*ER+pR)];

    % Add the matrix dissipation term to complete the Roe flux.
    Roe = ( FL + FR  - R*(ws.*dV))/2;
end

function Rusanov = RUSflux(qL,qR)
    % Rusanov flux 
    global gamma
    
    % Left state
    rL = qL(1);
    uL = qL(2)./qL(1);
    EL = qL(3)./rL;
    pL = (gamma-1)*( qL(3) - rL*uL*uL/2 );
    HL = ( qL(3) + pL )./ rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)./qR(1);
    ER = qR(3)./rR;
    pR = (gamma-1)*( qR(3) - rR*uR*uR/2 );
    HR = ( qR(3) + pR )./ rR;
    
    % First compute the Roe Averages
    RT = sqrt(rR/rL);
    %r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    H = (HL+RT*HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-u*u/2) );
    
    % Left and Right fluxes
    FL=[rL.*uL; rL.*uL.^2+pL; uL.*(rL.*EL+pL)];
    FR=[rR.*uR; rR.*uR.^2+pR; uR.*(rR.*ER+pR)];
    
    % Rusanov numerical flux
    smax = abs(u)+a;     Rusanov = ( FR + FL + smax*(qL-qR) )/2;
end

function AUSM = AUSMflux(qL,qR)
    % AUSM numerical flux
    %
    % M.-S. Liou and C. J. Steffen, A New Flux Splitting Scheme, Journal of
    % Computational Physics, 107, pp. 23-39, 1993.
    global gamma

    % Left state
    rL = qL(1);
    uL = qL(2)./qL(1);
    pL = (gamma-1)*( qL(3) - rL*uL*uL/2 );
    aL = sqrt(gamma*pL./rL);
    ML = uL/aL;
    HL = ( qL(3) + pL )./ rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)./qR(1);
    pR = (gamma-1)*( qR(3) - rR*uR*uR/2 );
    aR = sqrt(gamma*pR./rR);
    MR = uR/aR;
    HR = ( qR(3) + pR )./ rR;

    % Positive M and p in the LEFT cell.
    if (ML <= -1)
        Mp = 0;
        Pp = 0;
    elseif (ML < 1)
        Mp = (ML+1)*(ML+1)/4;
        Pp = pL*(1+ML)*(1+ML)*(2-ML)/4; % or use Pp = (1+ML)*pL/2
    else
        Mp = ML;
        Pp = pL;
    end

    % Negative M and p in the RIGHT cell.
    if   (MR <= -1)
        Mm = MR;
        Pm = pR;
    elseif (MR < 1)
        Mm = -(MR-1)*(MR-1)/4;
        Pm =  pR*(1-MR)*(1-MR)*(2+MR)/4; % or use Pm = (1-MR)*pR/2
    else
        Mm = 0;
        Pm = 0;
    end

    % Positive Part of Flux evaluated in the left cell.
    Fp(1) = max(0,Mp+Mm)*aL * rL;
    Fp(2) = max(0,Mp+Mm)*aL * rL*uL  + Pp;
    Fp(3) = max(0,Mp+Mm)*aL * rL*HL;

    % Negative Part of Flux evaluated in the right cell.
    Fm(1) = min(0,Mp+Mm)*aR * rR;
    Fm(2) = min(0,Mp+Mm)*aR * rR*uR  + Pm;
    Fm(3) = min(0,Mp+Mm)*aR * rR*HR;

    % Compute the flux: Fp(uL)+Fm(uR).
    AUSM = Fp + Fm;
end

function HLLE = HLLEflux(qL,qR)
    % Compute HLLE flux
    global gamma

    % Left state
    rL = qL(1);
    uL = qL(2)./rL;
    EL = qL(3)./rL;
    pL = (gamma-1)*( qL(3) - rL*uL*uL/2 );
    aL = sqrt(gamma*pL/rL);
    HL = ( qL(3) + pL )./ rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)./rR;
    ER = qR(3)./rR;
    pR = (gamma-1)*( qR(3) - rR*uR*uR/2 );
    aR = sqrt(gamma*pR/rR);
    HR = ( qR(3) + pR )./ rR;
    
    % Evaluate the two wave speeds: Einfeldt.
    RT = sqrt(rR/rL);
    u = (uL+RT*uR)/(1+RT);
    H = (HL+RT*HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-u*u/2) );
    
    % Wave speed estimates
    SLm = min(uL-aL, u-a);
    SRp = max(uR+aR, u+a);
    
    % Left and Right fluxes
    FL=[rL.*uL; rL.*uL.^2+pL; uL.*(rL.*EL+pL)];
    FR=[rR.*uR; rR.*uR.^2+pR; uR.*(rR.*ER+pR)];
    
    % Compute the HLL flux.
    if (0 <= SLm)  % Right-going supersonic flow
        HLLE = FL;
    elseif (SLm <= 0) && (0 <= SRp) % Subsonic flow
        select = 1;
        switch select
            case 1 % True HLLE function
                HLLE = ( SRp*FL - SLm*FR + SLm*SRp*(qR-qL) )/(SRp-SLm);
            case 2 % Rusanov flux ( as suggested by Toro's book )
                smax = max(abs(SLm),abs(SRp)); 
                HLLE = ( FR + FL + smax*(qL-qR) )/2; % Rusanov flux%
        end
    elseif  (0 >= SRp) % Left-going supersonic flow
        HLLE = FR;
    end
end

function HLLC = HLLCflux(qL,qR)
    % Compute HLLC flux
    global gamma

    % Left state
    rL = qL(1);
    uL = qL(2)./rL;
    EL = qL(3)./rL;
    pL = (gamma-1)*( qL(3) - rL*uL*uL/2 );
    aL = sqrt(gamma*pL/rL);
    
    % Right state
    rR = qR(1);
    uR = qR(2)./rR;
    ER = qR(3)./rR;
    pR = (gamma-1)*( qR(3) - rR*uR*uR/2 );
    aR = sqrt(gamma*pR/rR);
    
    % Left and Right fluxes
    FL=[rL.*uL; rL.*uL.^2+pL; uL.*(rL.*EL+pL)];
    FR=[rR.*uR; rR.*uR.^2+pR; uR.*(rR.*ER+pR)];

    % Compute guess pressure from PVRS Riemann solver
    PPV  = max(0 , 0.5*(pL+pR) + 0.5*(uL-uR) * (0.25*(rL+rR)*(aL+aR)));
    pmin = min(pL,pR);
    pmax = max(pL,pR);
    Qmax = pmax/pmin;
    Quser= 2.0; % <--- parameter manually set (I don't like this!)
    
     if (Qmax <= Quser) && (pmin <= PPV) && (PPV <= pmax)
     % Select PRVS Riemann solver
         pM = PPV;
         %uM = 0.5*(uL + uR) + 0.5*(pL - pR)/CUP;
      else
         if PPV < pmin
         % Select Two-Rarefaction Riemann solver
            PQ  = (pL/pR)^(gamma - 1.0)/(2.0*gamma);
            uM  = (PQ*uL/aL + uR/aR + 2/(gamma-1)*(PQ-1.0))/(PQ/aL+1.0/aR);
            PTL = 1 + (gamma-1)/2.0*(uL - uM)/aL;
            PTR = 1 + (gamma-1)/2.0*(uM - uR)/aR;
            pM  = 0.5*(pL*PTL^(2*gamma/(gamma-1)) + pR*PTR^(2*gamma/(gamma-1)));
         else 
         % Use Two-Shock Riemann solver with PVRS as estimate
            GEL = sqrt((2/(gamma+1)/rL)/((gamma-1)/(gamma+1)*pL + PPV));
            GER = sqrt((2/(gamma+1)/rR)/((gamma-1)/(gamma+1)*pR + PPV));
            pM  = (GEL*pL + GER*pR - (uR - uL))/(GEL + GER);
            %uM  = 0.5*(uL + uR) + 0.5*(GER*(pM - pR) - GEL*(pM - pL));
         end
      end

    % Estimate wave speeds: SL, SR and SM (Toro, 1994)
    if pM>pL; zL=sqrt(1+(gamma+1)/(2*gamma)*(pM/pL-1)); else, zL=1; end    
    if pM>pR; zR=sqrt(1+(gamma+1)/(2*gamma)*(pM/pR-1)); else, zR=1; end
  
	SL = uL - aL*zL;
    SR = uR + aR*zR;
    SM = (pL-pR + rR*uR*(SR-uR) - rL*uL*(SL-uL))/(rR*(SR-uR) - rL*(SL-uL));
    
    % Compute the HLL flux.
    if 0 <= SL  % Right-going supersonic flow
        HLLC = FL;
    elseif (SL <= 0) && (0 <= SM)	% Subsonic flow to the right
        qsL = rL*(SL-uL)/(SL-SM)*[1; SM; qL(3)/rL + (SM-uL)*(SM+pL/(rL*(SL-uL)))];
        HLLC = FL + SL*(qsL - qL);
    elseif (SM <= 0) && (0 <= SR)	% Subsonic flow to the Left
        qsR = rR*(SR-uR)/(SR-SM)*[1; SM; qR(3)/rR + (SM-uR)*(SM+pR/(rR*(SR-uR)))];
        HLLC = FR + SR*(qsR - qR);
    elseif  0 >= SR % Left-going supersonic flow
        HLLC = FR;
    end
end

function [wn,wp] = WENO5recon(w,N)
% *************************************************************************
% Input: u(i) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2)];
% Output: res = df/dx;
%
% Based on:
% C.W. Shu's Lectures notes on: 'ENO and WENO schemes for Hyperbolic
% Conservation Laws' 
%
% coded by Manuel Diaz, 02.10.2012, NTU Taiwan.
% last update on 2016.04.29, NHRI Taiwan.
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
%                           |___________S2__________|
%                           |                       |
%                   |___________S1__________|       |
%                   |                       |       |
%           |___________S0__________|       |       |
%         ..|---o---|---o---|---o---|---o---|---o---|...
%           | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                                  -|
%                                 i+1/2
%
%
%                   |___________S0__________|
%                   |                       |
%                   |       |___________S1__________|
%                   |       |                       |
%                   |       |       |___________S2__________|
%                 ..|---o---|---o---|---o---|---o---|---o---|...
%                   | I{i-1}|  I{i} | I{i+1}| I{i+2}| I{i+3}|
%                                   |+
%                                 i+1/2
%
% WENO stencil: S{i} = [ I{i-2},...,I{i+3} ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=3; I=R:(N-R); % R: stencil size

%% Right State Extrapolation $u_{i+1/2}^{-}$
vmm = w(:,I-2);
vm  = w(:,I-1);
v   = w(:, I );
vp  = w(:,I+1);
vpp = w(:,I+2);

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
wn  = w0n.*(2*vmm - 7*vm + 11*v)/6 ...
    + w1n.*( -vm  + 5*v  + 2*vp)/6 ...
    + w2n.*(2*v   + 5*vp - vpp )/6;

%% Left State Extrapolation $u_{i+1/2}^{+}$ 
umm = w(:,I-1);
um  = w(:, I );
u   = w(:,I+1);
up  = w(:,I+2);
upp = w(:,I+3);

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

% Numerical Flux at cell boundary, $u_{i+1/2}^{+}$;
wp  = w0p.*( -umm + 5*um + 2*u  )/6 ...
	+ w1p.*( 2*um + 5*u  - up   )/6 ...
	+ w2p.*(11*u  - 7*up + 2*upp)/6;
end

function [qn,qp] = WENO7recon(w,N)
% *************************************************************************
% Input: u(i) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2)];
% Output: res = df/dx;
%
% Based on:
% C.W. Shu's Lectures notes on: 'ENO and WENO schemes for Hyperbolic
% Conservation Laws' 
%
% coded by Manuel Diaz, 02.10.2012, NTU Taiwan.
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
%                               |_______________S3______________|
%                               |                               |
%                       |______________S2_______________|       |
%                       |                               |       |
%               |______________S1_______________|       |       |
%               |                               |       |       |
%       |_______________S0______________|       |       |       |
%     ..|---o---|---o---|---o---|---o---|---o---|---o---|---o---|...
%       | I{i-3}| I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}| I{i+3}|
%                                      -|
%                                     i+1/2
%
%       |______________S0_______________|
%       |                               |
%       |       |______________S1_______________|
%       |       |                               |
%       |       |       |______________S2_______________|
%       |       |       |                               |
%       |       |       |       |_______________S3______________|
%     ..|---o---|---o---|---o---|---o---|---o---|---o---|---o---|...
%       | I{i-3}| I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|| I{i+3}
%                               |+
%                             i-1/2
%
% WENO stencil: S{i} = [ I{i-3},...,I{i+3} ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The stencil size
I=4:(N-4); 

%% Right Flux: f_{i+1/2}^{-}
vmmm= w(:,I-3,:);
vmm = w(:,I-2,:);
vm  = w(:,I-1,:);
vo  = w(:, I ,:);
vp  = w(:,I+1,:);
vpp = w(:,I+2,:);
vppp= w(:,I+3,:);

% Smooth Indicators
B0n = vm.*(134241*vm-114894*vo)   +vmmm.*(56694*vm-47214*vmm+6649*vmmm-22778*vo)...
        +25729*vo.^2  +vmm.*(-210282*vm+85641*vmm+86214*vo);
B1n = vo.*(41001*vo-30414*vp)     +vmm.*(-19374*vm+3169*vmm+19014*vo-5978*vp)...
        +6649*vp.^2   +vm.*(33441*vm-70602*vo+23094*vp);
B2n = vp.*(33441*vp-19374*vpp)    +vm.*(6649*vm-30414*vo+23094*vp-5978*vpp)...
        +3169*vpp.^2  +vo.*(41001*vo-70602*vp+19014*vpp);
B3n = vpp.*(85641*vpp-47214*vppp) +vo.*(25729*vo-114894*vp+86214*vpp-22778*vppp)...
        +6649*vppp.^2 +vp.*(134241*vp-210282*vpp+56694*vppp);

% Constants
g0 = 1/35; g1 = 12/35; g2 = 18/35; g3 = 4/35; epsilon = 1e-6;

% Alpha weights
alpha0n = g0./(epsilon + B0n).^2;
alpha1n = g1./(epsilon + B1n).^2;
alpha2n = g2./(epsilon + B2n).^2;
alpha3n = g3./(epsilon + B3n).^2;
alphasumn = alpha0n + alpha1n + alpha2n + alpha3n;

% Non-linear weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
w2n = alpha2n./alphasumn;
w3n = alpha3n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
qn = w0n.*(-3*vmmm + 13*vmm - 23*vm  + 25*vo  )/12 + ...
     w1n.*( 1*vmm  - 5*vm   + 13*vo  +  3*vp  )/12 + ...
     w2n.*(-1*vm   + 7*vo   +  7*vp  -  1*vpp )/12 + ...
     w3n.*( 3*vo   + 13*vp  -  5*vpp +  1*vppp)/12;

%% Left Flux: f_{i+1/2}^{+}
ummm= w(:,I-2,:);
umm = w(:,I-1,:);
um  = w(:, I ,:);
uo  = w(:,I+1,:);
up  = w(:,I+2,:);
upp = w(:,I+3,:);
uppp= w(:,I+4,:);

% Smooth Indicators
B0p = um.*(134241*um-114894*uo)   +ummm.*(56694*um-47214*umm+6649*ummm-22778*uo)...
        +25729*uo.^2  +umm.*(-210282*um+85641*umm+86214*uo);
B1p = uo.*(41001*uo-30414*up)     +umm.*(-19374*um+3169*umm+19014*uo-5978*up)...
        +6649*up.^2   +um.*(33441*um-70602*uo+23094*up);
B2p = up.*(33441*up-19374*upp)    +um.*(6649*um-30414*uo+23094*up-5978*upp)...
        +3169*upp.^2  +uo.*(41001*uo-70602*up+19014*upp);
B3p = upp.*(85641*upp-47214*uppp) +uo.*(25729*uo-114894*up+86214*upp-22778*uppp)...
        +6649*uppp.^2 +up.*(134241*up-210282*upp+56694*uppp);

% Constants
g0 = 4/35; g1 = 18/35; g2 = 12/35; g3 = 1/35; epsilon = 1e-6;

% Alpha weights
alpha0p = g0./(epsilon + B0p).^2;
alpha1p = g1./(epsilon + B1p).^2;
alpha2p = g2./(epsilon + B2p).^2;
alpha3p = g3./(epsilon + B3p).^2;
alphasump = alpha0p + alpha1p + alpha2p + alpha3p;

% Non-linear weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;
w2p = alpha2p./alphasump;
w3p = alpha3p./alphasump;

% Numerical Flux at cell boundary, $u_{i+1/2}^{+}$;
qp = w0p.*( 1*ummm - 5*umm  + 13*um  +  3*uo  )/12 + ...
     w1p.*(-1*umm  + 7*um   +  7*uo  -  1*up  )/12 + ... 
     w2p.*( 3*um   + 13*uo  -  5*up  +  1*upp )/12 + ...
     w3p.*(25*uo   - 23*up  + 13*upp -  3*uppp)/12;
end