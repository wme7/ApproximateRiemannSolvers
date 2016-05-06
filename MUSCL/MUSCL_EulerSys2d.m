function [res] = MUSCL_EulerSys2d(q,smax,gamma,dx,dy,N,M,limiter,fluxMethod)
%   MUSCL Monotonic Upstreat Centered Scheme for Conservation Laws
%   Van Leer's MUSCL reconstruction scheme using piece wise linear
%   reconstruction
%  
%   e.g. where: limiter='MC'; fluxMethod='AUSM';
%
%   Flux at j+1/2
% 
%     j+1/2         Cell's grid:
%   | wL|   |
%   |  /|wR |           1   2   3   4        N-2 N-1  N
%   | / |\  |   {x=0} |-o-|-o-|-o-|-o-| ... |-o-|-o-|-o-| {x=L}
%   |/  | \ |         1   2   3   4   5        N-1  N  N+1
%   |   |  \|
%   |   |   |       NC: Here cells 1 and N are ghost cells
%     j  j+1            faces 2 and N, are the real boundary faces.
%
%   q = cat(3, r, r.*u, r.*v, r.*E);
%   F = cat(3, r.*u, r.*u.^2+p, r.*u.*v, u.*(r.*E+p));
%   G = cat(3, r.*v, r.*u.*v, r.*v.^2+p, v.*(r.*E+p));
%
% Written by Manuel Diaz, NTU, 04.29.2015.
    res = zeros(M,N,4);

    % Normal unitary face vectors: (nx,ny)
    %normals = {[0,1], [1,0], [0,-1], [-1,0]}; % i.e.: [N, E, S, W] 

    % Build cells
    cell(M,N).all = M*N;
    for i = 1:M
        for j = 1:N
            cell(i,j).q = [q(i,j,1);q(i,j,2);q(i,j,3);q(i,j,4)];
            cell(i,j).dqdy = zeros(4,1);
            cell(i,j).dqdx = zeros(4,1);
            cell(i,j).res = zeros(4,1);
        end
    end

    % Compute and limit slopes at cells (i,j)
    for i = 2:M-1       % only internal cells
        for j = 2:N-1   % only internal cells
            for k = 1:4;
                switch limiter
                    case 'MC' % MC limiter
                        % Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}
                        dqw = 2*(cell(i, j ).q(k) - cell(i,j-1).q(k))/dx; % du btn j and j-1
                        dqe = 2*(cell(i,j+1).q(k) - cell(i, j ).q(k))/dx; % du btn j+1 and j
                        dqc = (cell(i,j+1).q(k)-cell(i,j-1).q(k))/(2*dx); % du btn j+1 and j-1
                        cell(i,j).dqdx(k) = minmod([dqw,dqe,dqc]);
                        dqs = 2*(cell( i ,j).q(k) - cell(i-1,j).q(k))/dy; % du btn i and i-1
                        dqn = 2*(cell(i+1,j).q(k) - cell( i ,j).q(k))/dy; % du btn i+1 and i
                        dqc = (cell(i+1,j).q(k)-cell(i-1,j).q(k))/(2*dy); % du btn j+1 and j-1
                        cell(i,j).dqdy(k) = minmod([dqs,dqn,dqc]);
                    case 'MM' % Minmod limiter
                        % Find dq_j = minmod{fwd diff, bwd diff}
                        dqw = (cell(i, j ).q(k) - cell(i,j-1).q(k))/dx; % du btn j and j-1
                        dqe = (cell(i,j+1).q(k) - cell(i, j ).q(k))/dx; % du btn j+1 and j
                        cell(i,j).dqdx(k) = minmod([dqw,dqe]);
                        dqs = (cell( i ,j).q(k) - cell(i-1,j).q(k))/dy; % du btn i and i-1
                        dqn = (cell(i+1,j).q(k) - cell( i ,j).q(k))/dy; % du btn i+1 and i
                        cell(i,j).dqdy(k) = minmod([dqs,dqn]);
                    case 'VA' % Van Albada limiter
                        % Find dq_j = minmod{fwd diff, bwd diff}
                        dqw = (cell(i, j ).q(k) - cell(i,j-1).q(k))/dx; % du btn j and j-1
                        dqe = (cell(i,j+1).q(k) - cell(i, j ).q(k))/dx; % du btn j+1 and j
                        cell(i,j).dqdx(k) = vanalbada(dqw,dqe,dx);
                        dqs = (cell( i ,j).q(k) - cell(i-1,j).q(k))/dy; % du btn i and i-1
                        dqn = (cell(i+1,j).q(k) - cell( i ,j).q(k))/dy; % du btn i+1 and i
                        cell(i,j).dqdy(k) = vanalbada(dqs,dqn,dy);
                end
            end
        end
    end

	%%%%%%%%%%%%%
    % Residuals %
    %%%%%%%%%%%%%
    
    % Compute residuals x-direction
    for i = 2:M-1
        for j = 2:N-2
            % Left (inside) and Right (outside) extrapolated q-values at the boundaries
            qxL = cell( i,j ).q + cell( i,j ).dqdx*dx/2; % q_{i,j+1/2}^{-} from (i,j)
            qxR = cell(i,j+1).q - cell(i,j+1).dqdx*dx/2; % q_{i,j+1/2}^{+} from (i,j+1)
            % compute flux at j+1/2 using
            switch fluxMethod
                case 'LF'  % Lax-Friedrichs
                    flux = LFflux(qxL,qxR,gamma,[1,0],smax); % F_{i,j+1/2}
                case 'RUS' % Rusanov or local Lax-Friedrichs
                    flux = RUSflux(qxL,qxR,gamma,[1,0]);    % F_{i,j+1/2}
                case 'ROE' % Roe flux
                    flux = ROEflux(qxL,qxR,gamma,[1,0]);    % F_{i,j+1/2}
                case 'HLLE' % HLLE flux
                    flux = HLLEflux(qxL,qxR,gamma,[1,0]);   % F_{i,j+1/2}
               	case 'RHLLE' % RHLLE flux
                    flux = RHLLEflux(qxL,qxR,gamma,[1,0]);  % F_{i,j+1/2}
                otherwise
                    error('flux option not available')
            end
            % contributions to the residual of cell (i,j) and cells around it
            cell( i,j ).res = cell( i,j ).res + flux/dx;
            cell(i,j+1).res = cell(i,j+1).res - flux/dx;
        end
    end
    
    % Compute residuals y-direction
    for i = 2:M-2
        for j = 2:N-1
            % Left (inside) and Right (outside) extrapolated q-values at the boundaries
            qyL = cell( i,j ).q + cell( i,j ).dqdy*dy/2; % q_{i+1/2,j}^{-} from (i,j)
            qyR = cell(i+1,j).q - cell(i+1,j).dqdy*dy/2; % q_{i+1/2,j}^{+} from (i,j+1)
            % compute flux at j+1/2 using
            switch fluxMethod
                case 'LF'  % Lax-Friedrichs
                    flux = LFflux(qyL,qyR,gamma,[0,1],smax); % F_{i+1/2,j}
                case 'RUS' % Rusanov or local Lax-Friedrichs
                    flux = RUSflux(qyL,qyR,gamma,[0,1]);     % F_{i+1/2,j}
                case 'ROE' % Roe flux
                    flux = ROEflux(qyL,qyR,gamma,[0,1]);     % F_{i+1/2,j}
                case 'HLLE' % HLLE flux
                    flux = HLLEflux(qyL,qyR,gamma,[0,1]);    % F_{i+1/2,j}
                case 'RHLLE' % RHLLE flux
                    flux = RHLLEflux(qyL,qyR,gamma,[0,1]);   % F_{i+1/2,j}
                otherwise
                    error('flux option not available')
            end
            % contributions to the residual of cell (i,j) and cells around it
            cell( i,j ).res = cell( i,j ).res + flux/dy;
            cell(i+1,j).res = cell(i+1,j).res - flux/dy;
        end
    end
    
    %%%%%%%%%%%
    % set BCs %
    %%%%%%%%%%%
    
    % Flux contribution of the MOST NORTH FACE: north face of cells j=M-1.
    for j = 2:N-1
        qR = cell(M-1,j).q + cell(M-1,j).dqdy*dy/2;     qL = qR;
        switch fluxMethod
            case 'LF'  % Lax-Friedrichs
                flux = LFflux(qL,qR,gamma,[0,1],smax);	% F_{i+1/2,j}
            case 'RUS' % Rusanov or local Lax-Friedrichs
                flux = RUSflux(qL,qR,gamma,[0,1]);      % F_{i+1/2,j}
            case 'ROE' % Roe flux
                flux = ROEflux(qL,qR,gamma,[0,1]);      % F_{i+1/2,j}
            case 'HLLE' % HLLE flux
                flux = HLLEflux(qL,qR,gamma,[0,1]);     % F_{i+1/2,j}
          	case 'RHLLE' % RHLLE flux
                flux = RHLLEflux(qL,qR,gamma,[0,1]);    % F_{i+1/2,j}
        end
        cell(M-1,j).res = cell(M-1,j).res + flux/dy;
    end
    
    % Flux contribution of the MOST EAST FACE: east face of cell j=N-1.
    for i = 2:M-1
        qR = cell(i,N-1).q + cell(i,N-1).dqdx*dx/2;     qL = qR;
        switch fluxMethod
            case 'LF'  % Lax-Friedrichs
                flux = LFflux(qL,qR,gamma,[1,0],smax);	% F_{j,i+1/2}
            case 'RUS' % Rusanov or local Lax-Friedrichs
                flux = RUSflux(qL,qR,gamma,[1,0]);      % F_{i,j+1/2}
            case 'ROE' % Roe flux
                flux = ROEflux(qL,qR,gamma,[1,0]);      % F_{i,j+1/2}
            case 'HLLE' % HLLE flux
                flux = HLLEflux(qL,qR,gamma,[1,0]);     % F_{i,j+1/2}
          	case 'RHLLE' % RHLLE flux
                flux = RHLLEflux(qL,qR,gamma,[1,0]);    % F_{i,j+1/2}
        end
        cell(i,N-1).res = cell(i,N-1).res + flux/dx;
    end
    
    % Flux contribution of the MOST SOUTH FACE: south face of cells j=2.
    for j = 2:N-1
        qR = cell(2,j).q - cell(2,j).dqdy*dy/2;     qL = qR;
        switch fluxMethod
            case 'LF'  % Lax-Friedrichs
                flux = LFflux(qL,qR,gamma,[0,-1],smax); % F_{i-1/2,j}
            case 'RUS' % Rusanov or local Lax-Friedrichs
                flux = RUSflux(qL,qR,gamma,[0,-1]);     % F_{i-1/2,j}
            case 'ROE' % Roe flux
                flux = ROEflux(qL,qR,gamma,[0,-1]);     % F_{i-1/2,j}
            case 'HLLE' % HLLE flux
                flux = HLLEflux(qL,qR,gamma,[0,-1]);    % F_{i-1/2,j}
           	case 'RHLLE' % RHLLE flux
                flux = RHLLEflux(qL,qR,gamma,[0,-1]);   % F_{i-1/2,j}
        end
        cell(2,j).res = cell(2,j).res + flux/dy;
    end
    
    % Flux contribution of the MOST WEST FACE: west face of cells j=2.
    for i = 2:M-1
        qR = cell(i,2).q - cell(i,2).dqdx*dx/2;     qL = qR;
        switch fluxMethod
            case 'LF'  % Lax-Friedrichs
                flux = LFflux(qL,qR,gamma,[-1,0],smax); % F_{i,j-1/2}
            case 'RUS' % Rusanov or local Lax-Friedrichs
                flux = RUSflux(qL,qR,gamma,[-1,0]);     % F_{i,j-1/2}
            case 'ROE' % Roe flux
                flux = ROEflux(qL,qR,gamma,[-1,0]);     % F_{i,j-1/2}
            case 'HLLE' % HLLE flux
                flux = HLLEflux(qL,qR,gamma,[-1,0]);	% F_{i,j-1/2}
          	case 'RHLLE' % RHLLE flux
                flux = RHLLEflux(qL,qR,gamma,[-1,0]);	% F_{i,j-1/2}
        end
        cell(i,2).res = cell(i,2).res + flux/dx;
    end
    
    % Prepare residual as layers: [rho, rho*u, rho*v, rho*E]
    parfor i = 1:M
        for j = 1:N
            res(i,j,:) = cell(i,j).res;
        end
    end
end

function mm = minmod(v)
    % Using Harten's generalized definition
    % minmod: zero if opposite sign, otherwise the one of smaller magnitude.
    s = sum(sign(v))/numel(v); 
    if abs(s)==1; mm = s*min(abs(v(:))); else mm=0; end
    %m=size(v,1); mm=zeros(size(v,1),1); s=sum(sign(v),2)/m; ids=find(abs(s)==1);
    %if(~isempty(ids)); mm(ids)=s(ids).*min(abs(v(ids,:)),[],2); end
end

function va = vanalbada(da,db,h)
    % Van Albada Slope Limiter Function
    % vanAlbada: extend the simetric formulation of the van leer limiter
    eps2=(0.3*h)^3; 
    va=0.5*(sign(da)*sign(db)+1)*((db^2+eps2)*da+(da^2+eps2)*db)/(da^2+db^2+2*eps2);
end

function LF = LFflux(qL,qR,gamma,normal,smax)
    % Lax-Friedrichs flux

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

function Rusanov = RUSflux(qL,qR,gamma,normal)
    % Rusanov flux 
    
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

function Roe = ROEflux(qL,qR,gamma,normal)
    % Compute Roe flux
    
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

function HLLE = HLLEflux(qL,qR,gamma,normal)
    % Compute HLLE flux

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

function [RHLLE]= RHLLEflux(qL,qR,gamma,normal)
    % Evaluate Rotate-HLLE fluxes

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
    aL = sqrt(gamma*pL/rL);
    HL = ( qL(4) + pL ) / rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)/rR;
    vR = qR(3)/rR;
    vnR = uR*nx+vR*ny;
    vtR = uR*tx+vR*ty;
    pR = (gamma-1)*( qR(4) - rR*(uR^2+vR^2)/2 );
    aR = sqrt(gamma*pR/rR);
    HR = ( qR(4) + pR ) / rR;
    
    % Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
    % (NB: n1 and n2 may need to be frozen at some point during
    %     a steady calculation to fully make it converge. For time-accurate
    %     calculation, this is fine.)
    % NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).
    
    eps = 1.0e-12; % * M_inf;
    abs_dq = sqrt( (uR-uL)^2 + (vR-vL)^2 );
    
    if ( abs_dq > eps)
        nx1 = (uR-uL)/abs_dq;
        ny1 = (vR-vL)/abs_dq;
    else
        nx1 = -ny;
        ny1 =  nx;
    end
    
    % Rey = 1000.0_p2
    % temp = ( tanh(Rey*(abs_dq-eps)) - tanh(-Rey*eps) ) &
    %       /( tanh(Rey*(   one-eps)) - tanh(-Rey*eps) )
    % nx1 = temp*(uR-uL)/(abs_dq + eps) + (one-temp)*(-ny)
    % ny1 = temp*(vR-vL)/(abs_dq + eps) + (one-temp)*( nx)
    
    alpha1 = nx * nx1 + ny * ny1;
    % To make alpha1 always positive.
    temp = sign(alpha1);
    nx1 = temp * nx1;
    ny1 = temp * ny1;
    alpha1 = temp * alpha1;
    
    % Take n2 as perpendicular to n1.
    nx2 = -ny1;
    ny2 =  nx1;
    alpha2 = nx * nx2 + ny * ny2;
    % To make alpha2 always positive.
    temp = sign(alpha2);
    nx2 = temp * nx2;
    ny2 = temp * ny2;
    alpha2 = temp * alpha2;
    
    %Now we are going to compute the Roe flux with n2 as the normal
    %and n1 as the tagent vector, with modified wave speeds (5.12)
    
    % First compute the Roe Averages
    RT = sqrt(rR/rL);
    r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = ( HL+RT* HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-(u^2+v^2)/2) );
    vn = u*nx2+v*ny2;
    vt = u*nx1+v*ny1;
    
    % Wave Strengths
    dr = rR - rL;     dp = pR - pL;     dvn= vnR - vnL;     dvt= vtR - vtL;
    dV = [(dp-r*a*dvn )/(2*a^2); r*dvt/a; dr-dp/(a^2); (dp+r*a*dvn)/(2*a^2)];
    
    % Wave Speed
    ws = [vn-a; vn; vn; vn+a]; abs_ws = abs(ws);
    
    % Harten's Entropy Fix JCP(1983), 49, pp357-393:
    % only for the nonlinear fields.
    dws(1)=1/5; if ws(1)<dws(1); abs_ws(1)=( abs_ws(1)*abs_ws(1)/dws(1)+dws(1) )/2; end
    dws(4)=1/5; if ws(4)<dws(4); abs_ws(4)=( abs_ws(4)*abs_ws(4)/dws(4)+dws(4) )/2; end

    % Wave speed estimates, evaluated with [nx1,ny1] (= tangent wrt n2)
    SLm = min([ vtL-aL, vt-a, 0]);
    SRp = max([ vtR+aR, vt+a, 0]);
    
    % Modifed wave speed for the Rotated-RHLL flux (5.12) in the paper.
    ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + 2*alpha1*SRp*SLm )/ (SRp-SLm);

    %Right Eigenvectors: with n2 as normal and n1 as tangent.
    tx = nx1;
    ty = ny1;
    
    % Right Eigenvectors       
    Rv = [  1   ,  0   ,    1      ,  1   ;
          u-a*nx, a*tx ,    u      ,u+a*nx;
          u-a*ny, a*ty ,    u      ,u+a*ny;
          H-vn*a, vt*a ,(u^2+v^2)/2,H+vn*a];
    
    % Left and Right fluxes
    FL=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*ny; rL*vnL*HL];
    FR=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*ny; rR*vnR*HR];
    
    % Compute the HLL flux.
    RHLLE = ( SRp*FL - SLm*FR )/(SRp-SLm) - 0.5*Rv*(ws.*dV);
end