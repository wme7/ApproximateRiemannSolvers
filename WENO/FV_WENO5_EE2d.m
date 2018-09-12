function res = FV_WENO5_EE2d(q,smax,nx,ny,dx,dy,fluxMethod)
global gamma

% Compute primitive variables at solution points
w(:,:,1) = q(:,:,1);
w(:,:,2) = q(:,:,2)./q(:,:,1);
w(:,:,3) = q(:,:,3)./q(:,:,1);
w(:,:,4) = (gamma-1)*( q(:,:,4) - 0.5*(q(:,:,2).^2+q(:,:,3).^2)./q(:,:,1));

% 1. Reconstruct interface values: qL=q_{i+1/2}^{-} and qR=q_{i-1/2}^{+}
Recon = 'WENO5';
switch Recon
    case 'WENO5', [wL,wR] = WENO5recon_X(w,nx); R=3; % R: stencil size
    case 'WENO7', [wL,wR] = WENO7recon_X(w,nx); R=4;
    case 'Poly5', [wL,wR] = POLY5recon_X(w,nx); R=3;
    case 'Poly7', [wL,wR] = POLY7recon_X(w,nx); R=4;
    otherwise, error('reconstruction not available ;P');
end

% Compute conservative variables at faces
qL(:,:,1) = wL(:,:,1);
qL(:,:,2) = wL(:,:,2).*wL(:,:,1);
qL(:,:,3) = wL(:,:,3).*wL(:,:,1);
qL(:,:,4) = wL(:,:,4)./(gamma-1) + 0.5*wL(:,:,1).*(wL(:,:,2).^2+wL(:,:,3).^2);

qR(:,:,1) = wR(:,:,1);
qR(:,:,2) = wR(:,:,2).*wR(:,:,1);
qR(:,:,3) = wR(:,:,3).*wR(:,:,1);
qR(:,:,4) = wR(:,:,4)./(gamma-1) + 0.5*wR(:,:,1).*(wR(:,:,2).^2+wR(:,:,3).^2);

% 2. Compute finite volume residual term, df/dx.
res=zeros(size(q)); flux=zeros(size(q));

% Normal unitary face vectors: (nx,ny)
% normals = {[1,0], [0,1]}; % i.e.: x-axis, y-axis

for j=R:(ny-R)
    for i=R:(nx-R), I=i-2; % for all interior faces of the domain
        % compute flux at (i+1/2,j)
        switch fluxMethod
            case 'LF',  flux(j,i,:) = LFflux(squeeze(qL(j,I,:)),squeeze(qR(j,I,:)),[1,0],smax); % Lax Friedrichs
            case 'ROE', flux(j,i,:) = ROEflux(squeeze(qL(j,I,:)),squeeze(qR(j,I,:)),[1,0]); % Roe
            case 'RUS', flux(j,i,:) = RUSflux(squeeze(qL(j,I,:)),squeeze(qR(j,I,:)),[1,0]); % Rusanov
            case 'HLLE',flux(j,i,:) = HLLEflux(squeeze(qL(j,I,:)),squeeze(qR(j,I,:)),[1,0]); % HLLE
            case 'HLLC',flux(j,i,:) = HLLCflux(squeeze(qL(j,I,:)),squeeze(qR(j,I,:)),[1,0]); % HLLC
        end
        % Flux contribution to the residual of every cell
        res(j, i ,:) = res(j, i ,:) + flux(j,i,:)/dx;
        res(j,i+1,:) = res(j,i+1,:) - flux(j,i,:)/dx;
    end
end

% Flux contribution of the MOST WEST FACE: left face of cell j=3.
for j=R:(ny-R)
    switch fluxMethod
        case 'LF',  flux(j,3,:) = LFflux(squeeze(qL(j,R,:)),squeeze(qR(j,R,:)),[1,0],smax); % Lax Friedrichs
        case 'ROE', flux(j,3,:) = ROEflux(squeeze(qL(j,R,:)),squeeze(qR(j,R,:)),[1,0]); % Roe
        case 'RUS', flux(j,3,:) = RUSflux(squeeze(qL(j,R,:)),squeeze(qR(j,R,:)),[1,0]); % Rusanov
        case 'HLLE',flux(j,3,:) = HLLEflux(squeeze(qL(j,R,:)),squeeze(qR(j,R,:)),[1,0]); % HLLE
        case 'HLLC',flux(j,3,:) = HLLCflux(squeeze(qL(j,R,:)),squeeze(qR(j,R,:)),[1,0]); % HLLE
    end
    res(j,3,:) = res(j,3,:) - flux(j,3,:)/dx;
end

 
% Flux contribution of the MOST EAST FACE: right face of cell j=nx-2.
I=nx+1-2*R;
for j=R:(ny-R)
    switch fluxMethod
        case 'LF',  flux(j,nx-2,:) = LFflux(squeeze(qL(j,I,:)),squeeze(qR(j,I,:)),[1,0],smax); % Lax Friedrichs
        case 'ROE', flux(j,nx-2,:) = ROEflux(squeeze(qL(j,I,:)),squeeze(qR(j,I,:)),[1,0]); % Roe
        case 'RUS', flux(j,nx-2,:) = RUSflux(squeeze(qL(j,I,:)),squeeze(qR(j,I,:)),[1,0]); % Rusanov
        case 'HLLE',flux(j,nx-2,:) = HLLEflux(squeeze(qL(j,I,:)),squeeze(qR(j,I,:)),[1,0]); % HLLE
        case 'HLLC',flux(j,nx-2,:) = HLLCflux(squeeze(qL(j,I,:)),squeeze(qR(j,I,:)),[1,0]); % HLLC
    end
    res(j,nx-2,:) = res(j,nx-2,:) + flux(j,nx-2,:)/dx;
end

% 1. Reconstruct interface values: qL=q_{j+1/2}^{-} and qR=q_{j-1/2}^{+}
switch Recon
    case 'WENO5', [qL,qR] = WENO5recon_Y(q,ny); R=3; % R: stencil size
    case 'WENO7', [qL,qR] = WENO7recon_Y(q,ny); R=4;
    case 'Poly5', [qL,qR] = POLY5recon_Y(q,ny); R=3;
    case 'Poly7', [qL,qR] = POLY7recon_Y(q,ny); R=4;
    otherwise, error('reconstruction not available ;P');
end

%res=zeros(size(w)); flux=zeros(size(w));
for i=R:(nx-R)
    for j=R:(ny-R), J=j-2; % for all interior faces of the domain
        % compute flux at (i,j+1/2)
        switch fluxMethod
            case 'LF',  flux(j,i,:) = LFflux(squeeze(qL(J,i,:)),squeeze(qR(J,i,:)),[0,1],smax); % Lax Friedrichs
            case 'ROE', flux(j,i,:) = ROEflux(squeeze(qL(J,i,:)),squeeze(qR(J,i,:)),[0,1]); % Roe
            case 'RUS', flux(j,i,:) = RUSflux(squeeze(qL(J,i,:)),squeeze(qR(J,i,:)),[0,1]); % Rusanov
            case 'HLLE',flux(j,i,:) = HLLEflux(squeeze(qL(J,i,:)),squeeze(qR(J,i,:)),[0,1]); % HLLE
            case 'HLLC',flux(j,i,:) = HLLCflux(squeeze(qL(J,i,:)),squeeze(qR(J,i,:)),[0,1]); % HLLC
        end
        % Flux contribution to the residual of every cell
        res( j ,i,:) = res( j ,i,:) + flux(j,i,:)/dy;
        res(j+1,i,:) = res(j+1,i,:) - flux(j,i,:)/dy;
    end
end

% Flux contribution of the MOST SOUTH FACE: left face of cell i=3.
for i=R:(nx-R)
    switch fluxMethod
        case 'LF',  flux(3,i,:) = LFflux(squeeze(qL(R,i,:)),squeeze(qR(R,i,:)),[0,1],smax); % Lax Friedrichs
        case 'ROE', flux(3,i,:) = ROEflux(squeeze(qL(R,i,:)),squeeze(qR(R,i,:)),[0,1]); % Roe
        case 'RUS', flux(3,i,:) = RUSflux(squeeze(qL(R,i,:)),squeeze(qR(R,i,:)),[0,1]); % Rusanov
        case 'HLLE',flux(3,i,:) = HLLEflux(squeeze(qL(R,i,:)),squeeze(qR(R,i,:)),[0,1]); % HLLE
        case 'HLLC',flux(3,i,:) = HLLCflux(squeeze(qL(R,i,:)),squeeze(qR(R,i,:)),[0,1]); % HLLC
    end
    res(R,i,:) = res(R,i,:) - flux(3,i,:)/dy;
end

 
% Flux contribution of the MOST NORTH FACE: right face of cell i=ny-2.
J=ny+1-2*R;
for i=R:(nx-R)
    switch fluxMethod
        case 'LF',  flux(ny-2,i,:) = LFflux(squeeze(qL(J,i,:)),squeeze(qR(J,i,:)),[0,1],smax); % Lax Friedrichs
        case 'ROE', flux(ny-2,i,:) = ROEflux(squeeze(qL(J,i,:)),squeeze(qR(J,i,:)),[0,1]); % Roe
        case 'RUS', flux(ny-2,i,:) = RUSflux(squeeze(qL(J,i,:)),squeeze(qR(J,i,:)),[0,1]); % Rusanov
        case 'HLLE',flux(ny-2,i,:) = HLLEflux(squeeze(qL(J,i,:)),squeeze(qR(J,i,:)),[0,1]); % HLLE
        case 'HLLC',flux(ny-2,i,:) = HLLCflux(squeeze(qL(J,i,:)),squeeze(qR(J,i,:)),[0,1]); % HLLC
    end
    res(ny-2,i,:) = res(ny-2,i,:) + flux(ny-2,i,:)/dy;
end

end % FVM WENO

%%%%%%%%%%%%%%%%%
% Flux functions
%%%%%%%%%%%%%%%%%

function LF = LFflux(qL,qR,normal,smax)
    % Lax-Friedrichs flux
    global gamma

    % Normal vectors
    nx = normal(1);
    ny = normal(2);
    
    % Left state
    rL = qL(1);
    uL = qL(2)/rL;
    vL = qL(3)/rL;
    vnL= uL*nx + vL*ny;
    pL = (gamma-1)*( qL(4) - 0.5*rL*(uL^2+vL^2) );
    HL = ( qL(4) + pL ) / rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)/rR;
    vR = qR(3)/rR;
    vnR= uR*nx + vR*ny;
    pR = (gamma-1)*( qR(4) - 0.5*rR*(uR^2+vR^2) );
    HR = ( qR(4) + pR ) / rR;
    
    % Left and Right fluxes
    FL=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*ny; rL*vnL*HL];
    FR=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*ny; rR*vnR*HR];
    
    % Rusanov numerical flux
    LF = 0.5*( FR + FL + smax*(qL-qR) );
end

function Rusanov = RUSflux(qL,qR,normal)
    % Rusanov flux 
    global gamma
    
    % Normal vectors
    nx = normal(1);
    ny = normal(2);
    
    % Left state
    rL = qL(1);
    uL = qL(2)/rL;
    vL = qL(3)/rL;
    vnL= uL*nx + vL*ny;
    pL = (gamma-1)*( qL(4) - 0.5*rL*(uL^2+vL^2) );
    HL = ( qL(4) + pL ) / rL;

    % Right state
    rR = qR(1);
    uR = qR(2)/rR;
    vR = qR(3)/rR;
    vnR= uR*nx + vR*ny;
    pR = (gamma-1)*( qR(4) - 0.5*rR*(uR^2+vR^2) );
    HR = ( qR(4) + pR ) / rR;
    
    % First compute the Roe Averages
    RT = sqrt(rR/rL);
    %r= RT*rL;
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = ( HL+RT* HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-(u^2+v^2)/2) );
    
    % Left and Right fluxes
    FL=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*ny; rL*vnL*HL];
    FR=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*ny; rR*vnR*HR];
    
    % Rusanov numerical flux
    smax = abs(sqrt(u^2+v^2))+a; Rusanov = 0.5*( FR + FL + smax*(qL-qR) );
end

function Roe = ROEflux(qL,qR,normal)
    % Compute Roe flux
    global gamma
    
    % normal vectors
    nx = normal(1);
    ny = normal(2);
    
    % Tangent vectors
    tx = -ny;
    ty = nx;
    
    % Left state
    rL = qL(1);
    uL = qL(2)/rL;
    vL = qL(3)/rL;
    vnL = uL*nx+vL*ny;
    vtL = uL*tx+vL*ty;
    pL = (gamma-1)*( qL(4) - rL*(uL^2+vL^2)/2 );
    %aL = sqrt(gamma*pL/rL);
    HL = ( qL(4) + pL ) / rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)/rR;
    vR = qR(3)/rR;
    vnR = uR*nx+vR*ny;
    vtR = uR*tx+vR*ty;
    pR = (gamma-1)*( qR(4) - rR*(uR^2+vR^2)/2 );
    %aR = sqrt(gamma*pR/rR);
    HR = ( qR(4) + pR ) / rR;
    
    % First compute the Roe Averages
    RT = sqrt(rR/rL);
    r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = ( HL+RT* HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-(u^2+v^2)/2) );
    vn = u*nx+v*ny;
    vt = u*tx+v*ty;
    
    % Wave Strengths
    dr = rR - rL;     dp = pR - pL;     dvn= vnR - vnL;     dvt= vtR - vtL;
    dV = [(dp-r*a*dvn )/(2*a^2); r*dvt/a; dr-dp/(a^2); (dp+r*a*dvn)/(2*a^2)];
    
    % Wave Speed
    ws = [abs(vn-a); abs(vn); abs(vn); abs(vn+a)];
    
    % Harten's Entropy Fix JCP(1983), 49, pp357-393:
    % only for the nonlinear fields.
    dws(1)=1/5; if ws(1)<dws(1); ws(1)=( ws(1)*ws(1)/dws(1)+dws(1) )/2; end
    dws(4)=1/5; if ws(4)<dws(4); ws(4)=( ws(4)*ws(4)/dws(4)+dws(4) )/2; end
    
    % Right Eigenvectors       
    Rv = [  1   ,  0   ,    1      ,  1   ;
          u-a*nx, a*tx ,    u      ,u+a*nx;
          u-a*ny, a*ty ,    u      ,u+a*ny;
          H-vn*a, vt*a ,(u^2+v^2)/2,H+vn*a];
    
    % Left and Right fluxes
    FL=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*ny; rL*vnL*HL];
    FR=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*ny; rR*vnR*HR];
    
    % Dissipation Term
    Roe = (FL + FR - Rv*(ws.*dV))/2;

end

function HLLE = HLLEflux(qL,qR,normal)
    % Compute HLLE flux
    global gamma

    % normal vectors
    nx = normal(1);
    ny = normal(2);
       
    % Left state
    rL = qL(1);
    uL = qL(2)/rL;
    vL = qL(3)/rL;
    vnL = uL*nx+vL*ny;
    pL = (gamma-1)*( qL(4) - rL*(uL^2+vL^2)/2 );
    aL = sqrt(gamma*pL/rL);
    HL = ( qL(4) + pL ) / rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)/rR;
    vR = qR(3)/rR;
    vnR = uR*nx+vR*ny;
    pR = (gamma-1)*( qR(4) - rR*(uR^2+vR^2)/2 );
    aR = sqrt(gamma*pR/rR);
    HR = ( qR(4) + pR ) / rR;
    
    % First compute the Roe Averages
    RT = sqrt(rR/rL); % r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = ( HL+RT* HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-(u^2+v^2)/2) );
    vn = u*nx+v*ny;
    
    % Wave speed estimates
    SLm = min([ vnL-aL, vn-a, 0]);
    SRp = max([ vnR+aR, vn+a, 0]);
    
    % Left and Right fluxes
    FL=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*ny; rL*vnL*HL];
    FR=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*ny; rR*vnR*HR];
    
    % Compute the HLL flux.
    HLLE = ( SRp*FL - SLm*FR + SLm*SRp*(qR-qL) )/(SRp-SLm);
end

function HLLC = HLLCflux(qL,qR,normal)
    % Compute HLLC flux
    global gamma

    % normal vectors
    nx = normal(1);
    ny = normal(2);
    
    % Left state
    rL = qL(1);
    uL = qL(2)/rL;
    vL = qL(3)/rL;
    vnL = uL*nx+vL*ny;
    pL = (gamma-1)*( qL(4) - rL*(uL^2+vL^2)/2 );
    aL = sqrt(gamma*pL/rL);
    HL = ( qL(4) + pL ) / rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)/rR;
    vR = qR(3)/rR;
    vnR = uR*nx+vR*ny;
    pR = (gamma-1)*( qR(4) - rR*(uR^2+vR^2)/2 );
    aR = sqrt(gamma*pR/rR);
    HR = ( qR(4) + pR ) / rR;
       
    % Left and Right fluxes
    FL=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*ny; rL*vnL*HL];
    FR=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*ny; rR*vnR*HR];

    % Compute guess pressure from PVRS Riemann solver
    PPV  = max(0 , 0.5*(pL+pR) + 0.5*(vnL-vnR) * (0.25*(rL+rR)*(aL+aR)));
    pmin = min(pL,pR);
    pmax = max(pL,pR);
    Qmax = pmax/pmin;
    Quser= 2.0; % <--- parameter manually set (I don't like this!)
    
     if (Qmax <= Quser) && (pmin <= PPV) && (PPV <= pmax)
     % Select PRVS Riemann solver
         pM = PPV;
      else
         if PPV < pmin
         % Select Two-Rarefaction Riemann solver
            PQ  = (pL/pR)^(gamma - 1.0)/(2.0*gamma);
            uM  = (PQ*vnL/aL + vnR/aR + 2/(gamma-1)*(PQ-1.0))/(PQ/aL+1.0/aR);
            PTL = 1 + (gamma-1)/2.0*(vnL - uM)/aL;
            PTR = 1 + (gamma-1)/2.0*(uM - vnR)/aR;
            pM  = 0.5*(pL*PTL^(2*gamma/(gamma-1)) + pR*PTR^(2*gamma/(gamma-1)));
         else 
         % Use Two-Shock Riemann solver with PVRS as estimate
            GEL = sqrt((2/(gamma+1)/rL)/((gamma-1)/(gamma+1)*pL + PPV));
            GER = sqrt((2/(gamma+1)/rR)/((gamma-1)/(gamma+1)*pR + PPV));
            pM  = (GEL*pL + GER*pR - (vnR - vnL))/(GEL + GER);
         end
      end

    % Estimate wave speeds: SL, SR and SM (Toro, 1994)
    if pM>pL; zL=sqrt(1+(gamma+1)/(2*gamma)*(pM/pL-1)); else, zL=1; end    
    if pM>pR; zR=sqrt(1+(gamma+1)/(2*gamma)*(pM/pR-1)); else, zR=1; end
  
	SL = vnL - aL*zL;
    SR = vnR + aR*zR;
    SM = (pL-pR + rR*vnR*(SR-vnR) - rL*vnL*(SL-vnL))/(rR*(SR-vnR) - rL*(SL-vnL));
    
    % Compute the HLL flux.
    if (0 <= SL)  % Right-going supersonic flow
        HLLC = FL;
    elseif (SL <= 0) && (0 <= SM)	% Subsonic flow to the right
        qsL = rL*(SL-vnL)/(SL-SM)*[1; SM*nx+uL*abs(ny); SM*ny+vL*abs(nx); qL(4)/rL + (SM-vnL)*(SM+pL/(rL*(SL-vnL)))];
        HLLC = FL + SL*(qsL - qL);
    elseif (SM <= 0) && (0 <= SR)	% Subsonic flow to the Left
        qsR = rR*(SR-vnR)/(SR-SM)*[1; SM*nx+uR*abs(ny); SM*ny+vR*abs(nx); qR(4)/rR + (SM-vnR)*(SM+pR/(rR*(SR-vnR)))];
        HLLC = FR + SR*(qsR - qR);
    elseif	(0 >= SR) % Left-going supersonic flow
        HLLC = FR;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WENO and Polynomail reconstructions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qn,qp] = WENO5recon_X(w,N)
% *************************************************************************
% Input: u(i) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2)];
% Output: res = df/dx;
%
% Based on:
% Shu, Chi-Wang. "High order weighted essentially nonoscillatory schemes
% for convection dominated problems." SIAM review 51.1 (2009): 82-126.  
%
% coded by Manuel Diaz, 2016.04.29, NHRI Taiwan.
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
I=3:(N-3); % R:3 stencil size

%% Left State Extrapolation $u_{i+1/2}^{-}$
vmm = w(:,I-2,:);
vm  = w(:,I-1,:);
v   = w(:, I ,:);
vp  = w(:,I+1,:);
vpp = w(:,I+2,:);

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
qn  = w0n.*(2*vmm - 7*vm + 11*v)/6 ...
    + w1n.*( -vm  + 5*v  + 2*vp)/6 ...
    + w2n.*(2*v   + 5*vp - vpp )/6;

%% Right State Extrapolation $u_{i+1/2}^{+}$ 
umm = w(:,I-1,:);
um  = w(:, I ,:);
u   = w(:,I+1,:);
up  = w(:,I+2,:);
upp = w(:,I+3,:);

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
qp  = w0p.*( -umm + 5*um + 2*u  )/6 ...
	+ w1p.*( 2*um + 5*u  - up   )/6 ...
	+ w2p.*(11*u  - 7*up + 2*upp)/6;
end

function [qn,qp] = WENO5recon_Y(w,N)
% *************************************************************************
% Input: u(i) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2)];
% Output: res = df/dx;
%
% Based on:
% Shu, Chi-Wang. "High order weighted essentially nonoscillatory schemes
% for convection dominated problems." SIAM review 51.1 (2009): 82-126.  
%
% coded by Manuel Diaz, 2016.04.29, NHRI Taiwan.
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
I=3:(N-3); % R: stencil size

%% Left State Extrapolation $u_{i+1/2}^{-}$
vmm = w(I-2,:,:);
vm  = w(I-1,:,:);
v   = w( I ,:,:);
vp  = w(I+1,:,:);
vpp = w(I+2,:,:);

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
qn  = w0n.*(2*vmm - 7*vm + 11*v)/6 ...
    + w1n.*( -vm  + 5*v  + 2*vp)/6 ...
    + w2n.*(2*v   + 5*vp - vpp )/6;

%% Right State Extrapolation $u_{i+1/2}^{+}$ 
umm = w(I-1,:,:);
um  = w( I ,:,:);
u   = w(I+1,:,:);
up  = w(I+2,:,:);
upp = w(I+3,:,:);

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
qp  = w0p.*( -umm + 5*um + 2*u  )/6 ...
	+ w1p.*( 2*um + 5*u  - up   )/6 ...
	+ w2p.*(11*u  - 7*up + 2*upp)/6;
end

function [qn,qp] = POLY5recon_X(w,N)
% Direct Polynomial reconstruction
I=3:(N-3); % R:3 stencil size

%% Left Flux: f_{i+1/2}^{-}
vmm = w(:,I-2,:);
vm  = w(:,I-1,:);
v   = w(:, I ,:);
vp  = w(:,I+1,:);
vpp = w(:,I+2,:);

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
qn = ( 2*vmm - 13*vm + 47*v + 27*vp - 3*vpp)/60;

%% Right Flux: f_{i+1/2}^{+}
umm = w(:,I-1,:);
um  = w(:, I ,:);
u   = w(:,I+1,:);
up  = w(:,I+2,:);
upp = w(:,I+3,:);

% Numerical Flux at cell boundary, $u_{i+1/2}^{+}$;
qp = (-3*umm + 27*um + 47*u - 13*up + 2*upp)/60;
end

function [qn,qp] = POLY5recon_Y(w,N)
% The stencil size
I=3:(N-3); 

%% Left Flux: f_{i+1/2}^{-}
vmm = w(I-2,:,:);
vm  = w(I-1,:,:);
v   = w( I ,:,:);
vp  = w(I+1,:,:);
vpp = w(I+2,:,:);

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
qn = ( 2*vmm - 13*vm + 47*v + 27*vp - 3*vpp)/60;

%% Right Flux: f_{i+1/2}^{+}
umm = w(I-1,:,:);
um  = w( I ,:,:);
u   = w(I+1,:,:);
up  = w(I+2,:,:);
upp = w(I+3,:,:);

% Numerical Flux at cell boundary, $u_{i+1/2}^{+}$;
qp = (-3*umm + 27*um + 47*u - 13*up + 2*upp)/60;
end

function [qn,qp] = POLY7recon_X(w,N)
% The stencil size
I=4:(N-4);  

%% Left Flux: f_{i+1/2}^{-}
vmmm= w(:,I-3,:);
vmm = w(:,I-2,:);
vm  = w(:,I-1,:);
vo  = w(:, I ,:);
vp  = w(:,I+1,:);
vpp = w(:,I+2,:);
vppp= w(:,I+3,:);

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
qn = (-3*vmmm + 25*vmm - 101*vm  + 319*vo + 214*vp - 38*vpp + 4*vppp)/420;

%% Right Flux: f_{i+1/2}^{+}
ummm= w(:,I-2,:);
umm = w(:,I-1,:);
um  = w(:, I ,:);
uo  = w(:,I+1,:);
up  = w(:,I+2,:);
upp = w(:,I+3,:);
uppp= w(:,I+4,:);

% Numerical Flux at cell boundary, $u_{i+1/2}^{+}$;
qp = (4*ummm - 38*umm  + 214*um  + 319*uo - 101*up + 25*upp - 3*uppp)/420;
end

function [qn,qp] = POLY7recon_Y(w,N)
% The stencil size
I=4:(N-4);  

%% Left Flux: f_{i+1/2}^{-}
vmmm= w(I-3,:,:);
vmm = w(I-2,:,:);
vm  = w(I-1,:,:);
vo  = w( I ,:,:);
vp  = w(I+1,:,:);
vpp = w(I+2,:,:);
vppp= w(I+3,:,:);

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
qn = (-3*vmmm + 25*vmm - 101*vm  + 319*vo + 214*vp - 38*vpp + 4*vppp)/420;

%% Right Flux: f_{i+1/2}^{+}
ummm= w(I-2,:,:);
umm = w(I-1,:,:);
um  = w( I ,:,:);
uo  = w(I+1,:,:);
up  = w(I+2,:,:);
upp = w(I+3,:,:);
uppp= w(I+4,:,:);

% Numerical Flux at cell boundary, $u_{i+1/2}^{+}$;
qp = (4*ummm - 38*umm  + 214*um  + 319*uo - 101*up + 25*upp - 3*uppp)/420;
end