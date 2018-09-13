function res = FV_WENO_EE2d_Cprototype(q,a,nx,ny,dx,dy,t,fluxMethod,Recon,Test)
% Compute RHS of the semi-discrete form of the Euler equations.
global gamma preshock postshock mesh_wedge_position

%   Flux at j+1/2
% 
%     j+1/2    Cell's grid: (assuming WENO5, R=3)
%   |   |   |                   {x=0}             {x=L}
%   | wL|   |                     |                 |
%   |  /|wR |           1   2   3 | 4   5        N-3|N-2 N-1  N
%   | / |\  |         |-o-|-o-|-o-|-o-|-o-| ... |-o-|-o-|-o-|-o-|---> j
%   |/  | \ |             1   2   3   4   6    N-4 N-3 N-2 N-1  
%   |   |  \|                    {1} {2} {3}  ...  {nf}
%   |   |   |       NC: Here cells 1 to 3 and N-2 to N are ghost cells
%     j  j+1            faces 3 and N-3, are the real boundary faces.
%
%   q = cat(3, r, ru, rv, E);
%   F = cat(3, ru, ru^2+p, ruv, ru(E+p));
%   G = cat(3, rv, ruv, rv^2+p, rv(E+p));

% Identify number of gost cells
switch Recon
    case {'WENO5','Poly5'}, R=3; % R: stencil size and number of gost cells
    otherwise, error('reconstruction not available ;P');
end

switch Test
    case 'Smooth' % Set Periodic BCs
        for i=1:R
            q(:,i,:)=q(:,nx-R+i,:); q(:,nx-2*R+i,:)=q(:,R+i,:);	% Periodic BCs
        end
        for j=1:R
            q(j,:,:)=q(ny-R+i,:,:); q(ny-2*R+j,:,:)=q(R+j,:,:);	% Periodic BCs
        end
    case 'Riemann' % Set outflow BCs
        for i=1:R
            q(:,i,:)=q(:,R+1,:); q(:,nx+1-i,:)=q(:,nx-R,:);	% Neumann BCs
        end
        for j=1:R
            q(j,:,:)=q(R+1,:,:); q(ny+1-j,:,:)=q(ny-R,:,:);	% Neumann BCs
        end
    case 'DMR' % Set DMR test BCs
        % Static BCs
        for j=1:ny
            for i=1:R
                q(j,i,:)=postshock; q(j,nx+1-i,:)=preshock; % Dirichlet BCs
            end
        end
        % Static BCs at the bottom of domain
        for i=1:R
            for j=R+1:nx-R
                if j<(mesh_wedge_position+R+1)
                    q(i,j,:)=postshock; % Dirichlet BCs
                else
                    q(i,j,:)=q(R+1,j,:); q(i,j,3)=-q(R+1,j,3); % EE reflective BC
                end
            end
        end
        % Time dependent BCs at the top of domain: moving shock
        for i=ny+1-R:ny % only gosht cells at the top
            for j=R+1:nx-R % evaluate all x domain
                if distance_to_shock(j*dx+dx/2,i*dy+dy/2,t) < 0 % mesh_shock
                    q(i,j,:)=postshock; % Dirichlet BCs
                else
                    q(i,j,:)=preshock; % Dirichlet BCs
                end
            end
        end
    otherwise, error('Test boundaries not set!');
end

% Compute primitive variables at solution points
w(:,:,1) = q(:,:,1);
w(:,:,2) = q(:,:,2)./q(:,:,1);
w(:,:,3) = q(:,:,3)./q(:,:,1);
w(:,:,4) = (gamma-1)*( q(:,:,4) - 0.5*(q(:,:,2).^2+q(:,:,3).^2)./q(:,:,1));

% 1. Reconstruct interface values: qL=q_{i+1/2}^{-} and qR=q_{i-1/2}^{+}
switch Recon
    case 'WENO5', [wL,wR] = WENO5recon(w(R+1:ny-R,:,:),nx);
    case 'Poly5', [wL,wR] = POLY5recon(w(R+1:ny-R,:,:),nx);
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
res=zeros(size(q)); flux=zeros(size(qR)); nc=ny-2*R; nf=nx+1-2*R;

% Normal unitary face vectors: (nx,ny)
% normals = {[1,0], [0,1]}; % i.e.: x-axis, y-axis

% compute flux at (i+1/2,j)
for j=1:nc % for all interior cells
    for i=1:nf % for all interior faces 
        switch fluxMethod
            case 'LF',  flux(j,i,:) = LFflux(squeeze(qL(j,i,:)),squeeze(qR(j,i,:)),[1,0],a); % Lax Friedrichs
            case 'ROE', flux(j,i,:) = ROEflux(squeeze(qL(j,i,:)),squeeze(qR(j,i,:)),[1,0]);  % Roe
            case 'RUS', flux(j,i,:) = RUSflux(squeeze(qL(j,i,:)),squeeze(qR(j,i,:)),[1,0]);  % Rusanov
            case 'HLLE',flux(j,i,:) = HLLEflux(squeeze(qL(j,i,:)),squeeze(qR(j,i,:)),[1,0]); % HLLE
            case 'HLLC',flux(j,i,:) = HLLCflux(squeeze(qL(j,i,:)),squeeze(qR(j,i,:)),[1,0]); % HLLC
            otherwise, error('flux method not available ;P');
        end
    end
end

% Flux contribution to the residual of every cell
for j=1:nc % for all interior cells
    res(j+R,R+1,:) = res(j+R,R+1,:) - flux(j,1,:)/dx; % left face of cell j=4.
    for i = 2:nf-1 % for all interior faces
        res(j+R,i+R-1,:) = res(j+R,i+R-1,:) + flux(j,i,:)/dx;
        res(j+R, i+R ,:) = res(j+R, i+R ,:) - flux(j,i,:)/dx;
    end
    res(j+R,nx-R,:) = res(j+R,nx-R,:) + flux(j,nf,:)/dx; % right face of cell j=N-3.
end

% Compute primitive variables at solution points (they still in memory)
clear qL qR flux;

% 1. Reconstruct interface values: qL=q_{j+1/2}^{-} and qR=q_{j-1/2}^{+}
switch Recon
    case 'WENO5', [wL,wR] = WENO5recon(w(:,R+1:nx-R,:),ny);
    case 'Poly5', [wL,wR] = POLY5recon(w(:,R+1:nx-R,:),ny);
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
flux=zeros(size(qL)); nc=nx-2*R; nf=ny+1-2*R;

% compute flux at (i,j+1/2)
for i=1:nc % for all interior cells
    for j=1:nf % for all interior faces 
        switch fluxMethod
            case 'LF',  flux(j,i,:) = LFflux(squeeze(qL(j,i,:)),squeeze(qR(j,i,:)),[0,1],a); % Lax Friedrichs
            case 'ROE', flux(j,i,:) = ROEflux(squeeze(qL(j,i,:)),squeeze(qR(j,i,:)),[0,1]);  % Roe
            case 'RUS', flux(j,i,:) = RUSflux(squeeze(qL(j,i,:)),squeeze(qR(j,i,:)),[0,1]);  % Rusanov
            case 'HLLE',flux(j,i,:) = HLLEflux(squeeze(qL(j,i,:)),squeeze(qR(j,i,:)),[0,1]); % HLLE
            case 'HLLC',flux(j,i,:) = HLLCflux(squeeze(qL(j,i,:)),squeeze(qR(j,i,:)),[0,1]); % HLLC
        end
    end
end

% Flux contribution to the residual of every cell
for i=1:nc % for all interior cells
    res(R+1,i+R,:) = res(R+1,i+R,:) - flux(1,i,:)/dy;
    for j=2:nf-1 % for all interior cells
        res(j+R-1,i+R,:) = res(j+R-1,i+R,:) + flux(j,i,:)/dy;
        res( j+R ,i+R,:) = res( j+R ,i+R,:) - flux(j,i,:)/dy;
    end
    res(ny-R,i+R,:) = res(ny-R,i+R,:) + flux(nf,i,:)/dy;
end

end % FVM WENO

%%%%%%%%%%%%%%%%%%%
% Distance to shock (for Double Mach Reflection)
%%%%%%%%%%%%%%%%%%%

function distance = distance_to_shock(x,y,t)
    global shock_speed
    shock_slope = 1/tan(pi/6); % from problem definition
    wedge_position = 1/6; % from problem definition
    distance = (shock_slope*(x-wedge_position-shock_speed*t)-y) / sqrt((shock_slope)^2+1);
end

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

function [qn,qp] = WENO5recon(w)
% WENO5 reconstruction
% *************************************************************************
% Input: u(i) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2) u(i+3)];
% Output: u(i+1/2)^{-} & u(i+1/2)^{+}
%
% coded by Manuel Diaz, 2016.04.29, NHRI Taiwan.
% *************************************************************************


end



function [qn,qp] = POLY5recon(w)
% Direct Polynomial reconstruction
% *************************************************************************
% Input: [u] = [u(i-2) u(i-1) u(i) u(i+1) u(i+2) u(i+3)];
% Output: u(i+1/2)^{-} & u(i+1/2)^{+}
%
% coded by Manuel Diaz, 2016.04.29, NHRI Taiwan.
% *************************************************************************
I=3 % R:3 stencil size

%% Left Flux: f_{i+1/2}^{-}
vmm = w(I-2);
vm  = w(I-1);
v   = w( I );
vp  = w(I+1);
vpp = w(I+2);

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
qn = ( 2*vmm - 13*vm + 47*v + 27*vp - 3*vpp)/60;

%% Right Flux: f_{i+1/2}^{+}
umm = w(I-1);
um  = w( I );
u   = w(I+1);
up  = w(I+2);
upp = w(I+3);

% Numerical Flux at cell boundary, $u_{i+1/2}^{+}$;
qp = (-3*umm + 27*um + 47*u - 13*up + 2*upp)/60;
end