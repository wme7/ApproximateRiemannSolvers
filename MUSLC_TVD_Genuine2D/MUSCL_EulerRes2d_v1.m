function [res] = MUSCL_EulerRes2d_v1(q,dt,dx,dy,N,M,limiter,fluxMethod)
%   A genuine 2d HLLE Riemnan solver for Euler Equations using a Monotonic
%   Upstreat Centered Scheme for Conservation Laws (MUSCL).
%  
%   e.g. where: limiter='MC'; fluxMethod='HLLE1d';
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
% Written by Manuel Diaz, NTU, 05.25.2015.
    res = zeros(M,N,4);

    % Normal unitary face vectors: (nx,ny)
    %normals = {[0,1], [1,0], [0,-1], [-1,0]}; % i.e.: [N, E, S, W] 
    
    % Build cells
    C(M,N).q = zeros(4,1);
    for i = 1:M
        for j = 1:N
            C(i,j).q = [q(i,j,1);q(i,j,2);q(i,j,3);q(i,j,4)];
            C(i,j).dqdy = zeros(4,1);
            C(i,j).dqdx = zeros(4,1);
            C(i,j).res = zeros(4,1);
        end
    end

    % Compute and limit slopes at cells (i,j)
    for i = 2:M-1       % only internal cells
        for j = 2:N-1   % only internal cells
            for k = 1:4
                dqw = 2*(C( i,j ).q(k) - C(i,j-1).q(k))/dx; % du btn j and j-1
                dqe = 2*(C(i,j+1).q(k) - C( i,j ).q(k))/dx; % du btn j+1 and j
                dqs = 2*(C( i,j ).q(k) - C(i-1,j).q(k))/dy; % du btn i and i-1
                dqn = 2*(C(i+1,j).q(k) - C( i,j ).q(k))/dy; % du btn i+1 and i
                switch limiter
                    case 'MC' % MC limiter
                        % Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}
                        dqc = (C(i,j+1).q(k)-C(i,j-1).q(k))/(2*dx); % du btn j+1 and j-1
                        C(i,j).dqdx(k) = minmod([dqw,dqe,dqc]);
                        dqc = (C(i+1,j).q(k)-C(i-1,j).q(k))/(2*dy); % du btn j+1 and j-1
                        C(i,j).dqdy(k) = minmod([dqs,dqn,dqc]);
                    case 'MM' % Minmod limiter
                        % Find dq_j = minmod{fwd diff, bwd diff}
                        C(i,j).dqdx(k) = minmod([dqw,dqe]);
                        C(i,j).dqdy(k) = minmod([dqs,dqn]);
                    case 'VA' % Van Albada limiter
                        % Find dq_j = vanAlvada{fwd diff, bwd diff, h }
                        C(i,j).dqdx(k) = vanalbada(dqw,dqe,dx);
                        C(i,j).dqdy(k) = vanalbada(dqs,dqn,dy);
                    case 'VL' % Van Leer limiter
                        % Find dq_j = vanAlvada{fwd diff, bwd diff, h }
                        C(i,j).dqdx(k) = vanLeer(dqw,dqe);
                        C(i,j).dqdy(k) = vanLeer(dqs,dqn);
                end
            end
        end
    end

	%%%%%%%%%%%%%
    % Residuals %
    %%%%%%%%%%%%%
    
    % Compute residuals 
    for i = 2:M-1
        for j = 2:N-2
            % Compute fluxes
            %[f2N,f2S,g2E,g2W,f1N,f1S,g1E,g1W,sN,sS,sE,sW]
            [~,D1f2S,~,~,~,D1f1S,~,~,D1Sn,D1Ss,~,~] = ... % Dual cell 1
                HLLE2d(C(i,j),C(i,j+1),C(i+1,j),C(i+1,j+1),gamma,dx,dy);
            
            [D2f2N,~,~,~,D2f1N,~,~,~,D2Sn,D2Ss,~,~] = ... % Dual cell 2
                HLLE2d(C(i-1,j),C(i-1,j+1),C(i,j),C(i,j+1),gamma,dx,dy);
            
            switch assembly
                case 'manual'   % Manually assemble the fluxes
                    thyD1 = dt/(2*dy)*max(D1Ss,D1Sn);
                    thyD2 = dt/(2*dy)*max(D2Ss,D2Sn);
                    thy = 1 - thyD1 - thyD2;
                    flux = thyD1*D1f2S + thy*(D1f1S+D2f1N)/2 + thyD2*D2f2N;
                case 'simpson'  % Using simpsons rule
                    flux = 0;
                otherwise
                    error('not a valid assemble :/');
            end
            
            % Contributions to the residual of cell (i,j) and cells around it
            C(i, j ).res = C( i , j ).res + flux/dx;
            C(i,j+1).res = C( i ,j+1).res - flux/dx;
        end
    end
    
    for i = 2:M-2
        for j = 2:N-1
            % Compute fluxes
            [~,~,D1g2E,~,~,~,D1g1E,~,~,~,D1Se,D1Sw] = ... % Dual cell 1
                HLLE2d(C(i,j),C(i,j+1),C(i+1,j),C(i+1,j+1),gamma,dx,dy);
            
            [~,~,~,D3g2W,~,~,~,D3g1W,~,~,D3Se,D3Sw] = ... % Dual cell 3
                HLLE2d(C(i,j-1),C(i,j),C(i+1,j),C(i+1,j-1),gamma,dx,dy);
            
            switch assembly
                case 'manual'   % Manually assemble the fluxes
                    thxD1 = dt/(2*dx)*max(D1Se,D1Sw);
                    thxD3 = dt/(2*dx)*max(D3Se,D3Sw);
                    thx = 1 - thxD1 - thxD3;
                    flux = thxD3*D3g2W + thx*(D1g1E+D3g1W)/2 + thxD1*D1g2E;
                case 'simpson'  % Using simpsons rule
                    flux = 0;
                otherwise
                    error('not a valid assemble :/');
            end
            
            % Contributions to the residual of cell (i,j) and cells around it
            C( i ,j).res = C( i ,j).res + flux/dy;
            C(i+1,j).res = C(i+1,j).res - flux/dy;
        end
    end
    
    %%%%%%%%%%%
    % set BCs %
    %%%%%%%%%%%
    
    % Flux contribution of the MOST NORTH FACE: north face of cells i=M-1.
    for j = 2:N-1
        % Compute fluxes
        [~,~,D1g2E,~,~,~,D1g1E,~,~,~,D1Se,D1Sw] = ... % Dual cell 1
            HLLE2d(C(M-1,j),C(M-1,j+1),C(M,j),C(M,j+1),gamma,dx,dy);

        [~,~,~,D3g2W,~,~,~,D3g1W,~,~,D3Se,D3Sw] = ... % Dual cell 3
            HLLE2d(C(M-1,j-1),C(M-1,j),C(M,j),C(M,j-1),gamma,dx,dy);

        switch assembly
            case 'manual'   % Manually assemble the fluxes
                thxD1 = dt/(2*dx)*max(D1Se,D1Sw);
                thxD3 = dt/(2*dx)*max(D3Se,D3Sw);
                thx = 1 - thxD1 - thxD3;
                flux = thxD3*D3g2W + thx*(D1g1E+D3g1W)/2 + thxD1*D1g2E;
            case 'simpson'  % Using simpsons rule
                flux = 0;
            otherwise
                error('not a valid assemble :/');
        end

        % Contributions to the residual of cell (i,j) and cells around it
        C(M-1,j).res = C(M-1,j).res + flux/dy;
    end
     
    % Flux contribution of the MOST SOUTH FACE: south face of cells i=2.
    for j = 2:N-1
        % Compute fluxes
        [~,~,D1g2E,~,~,~,D1g1E,~,~,~,D1Se,D1Sw] = ... % Dual cell 1
            HLLE2d(C(1,j),C(1,j+1),C(2,j),C(2,j+1),gamma,dx,dy);

        [~,~,~,D3g2W,~,~,~,D3g1W,~,~,D3Se,D3Sw] = ... % Dual cell 3
            HLLE2d(C(1,j-1),C(1,j),C(2,j),C(2,j-1),gamma,dx,dy);

        switch assembly
            case 'manual'   % Manually assemble the fluxes
                thxD1 = dt/(2*dx)*max(D1Se,D1Sw);
                thxD3 = dt/(2*dx)*max(D3Se,D3Sw);
                thx = 1 - thxD1 - thxD3;
                flux = thxD3*D3g2W + thx*(D1g1E+D3g1W)/2 + thxD1*D1g2E;
            case 'simpson'  % Using simpsons rule
                flux = 0;
            otherwise
                error('not a valid assemble :/');
        end

        % Contributions to the residual of cell (i,j) and cells around it
        C(2,j).res = C(2,j).res - flux/dy;
    end

    % Flux contribution of the MOST EAST FACE: east face of cell j=N-1.
    for i = 2:M-1
        [~,D1f2S,~,~,~,D1f1S,~,~,D1Sn,D1Ss,~,~] = ... % Dual cell 1
            HLLE2d(C(i,N-1),C(i,N),C(i+1,N-1),C(i+1,N),gamma,dx,dy);

        [D2f2N,~,~,~,D2f1N,~,~,~,D2Sn,D2Ss,~,~] = ... % Dual cell 2
            HLLE2d(C(i-1,N-1),C(i-1,N),C(i,N-1),C(i,N),gamma,dx,dy);

        switch assembly
            case 'manual'   % Manually assemble the fluxes
                thyD1 = dt/(2*dy)*max(D1Ss,D1Sn);
                thyD2 = dt/(2*dy)*max(D2Ss,D2Sn);
                thy = 1 - thyD1 - thyD2;
                flux = thyD1*D1f2S + thy*(D1f1S+D2f1N)/2 + thyD2*D2f2N;
            case 'simpson'  % Using simpsons rule
                flux = 0;
            otherwise
                error('not a valid assemble :/');
        end
        C(i,N-1).res = C(i,N-1).res + flux/dx;
    end

    % Flux contribution of the MOST WEST FACE: west face of cells j=2.
    for i = 2:M-1
        [~,D1f2S,~,~,~,D1f1S,~,~,D1Sn,D1Ss,~,~] = ... % Dual cell 1
            HLLE2d(C(i,1),C(i,2),C(i+1,1),C(i+1,2),gamma,dx,dy);

        [D2f2N,~,~,~,D2f1N,~,~,~,D2Sn,D2Ss,~,~] = ... % Dual cell 2
            HLLE2d(C(i-1,1),C(i-1,2),C(i,1),C(i,2),gamma,dx,dy);

        switch assembly
            case 'manual'   % Manually assemble the fluxes
                thyD1 = dt/(2*dy)*max(D1Ss,D1Sn);
                thyD2 = dt/(2*dy)*max(D2Ss,D2Sn);
                thy = 1 - thyD1 - thyD2;
                flux = thyD1*D1f2S + thy*(D1f1S+D2f1N)/2 + thyD2*D2f2N;
            case 'simpson'  % Using simpsons rule
                flux = 0;
            otherwise
                error('not a valid assemble :/');
        end
        C(i,2).res = C(i,2).res - flux/dx;
    end
    
    % Prepare residual as layers: [rho, rho*u, rho*v, rho*E]
    parfor i = 1:M
        for j = 1:N
            res(i,j,:) = C(i,j).res;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions %
%%%%%%%%%%%%%%%%%%%%%%%

function mm = minmod(v)
    % Using Harten's generalized definition
    % minmod: zero if opposite sign, otherwise the one of smaller magnitude.
    s=sum(sign(v))/numel(v); if abs(s)==1; mm=s*min(abs(v(:))); else, mm=0; end
    %m=size(v,1); mm=zeros(size(v,1),1); s=sum(sign(v),2)/m; ids=find(abs(s)==1);
    %if(~isempty(ids)); mm(ids)=s(ids).*min(abs(v(ids,:)),[],2); end
end

function va = vanalbada(da,db,h)
    % Van Albada Slope Limiter Function
    % vanAlbada: extend the simetric formulation of the van leer limiter
    eps2=(0.3*h)^3; 
    va=0.5*(sign(da)*sign(db)+1)*((db^2+eps2)*da+(da^2+eps2)*db)/(da^2+db^2+2*eps2);
end

function vl = vanLeer(da,db)
    % Van Leer Slope Limiter Function
    vl = 0; if bd~=0, r=da/db; vl=(r+abs(r))/(1+abs(r)); end
end

function [f,g] = Fluxes(q,gamma)
    % q state
    r = q(1);
    u = q(2)/r;
    v = q(3)/r;
    p = (gamma-1)*( q(4) - 0.5*r*(u^2+v^2) );
    H = ( q(4)+p )/r;
    
    % compute fluxes
    f=[r*u; r*u*u + p; r*u*v; r*u*H]; % f: x-flux
    g=[r*v; r*v*u; r*v*v + p; r*v*H]; % g: y-flux
end

function [qS,nF,tF,SLm,SRp] = HLLE1d(qL,qR,gamma,normal)
    % Compute 1-d HLLE normal and perpendicular fluxes

    % normal vectors
    nx = normal(1);
    ny = normal(2);
       
    % Left state
    rL = qL(1);
    uL = qL(2)/rL;
    vL = qL(3)/rL;
    vnL = uL*nx+vL*ny;
    pL = (gamma-1)*( qL(4) - 0.5*rL*(uL^2+vL^2) );
    aL = sqrt(gamma*pL/rL);
    HL = ( qL(4) + pL ) / rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)/rR;
    vR = qR(3)/rR;
    vnR = uR*nx+vR*ny;
    pR = (gamma-1)*( qR(4) - 0.5*rR*(uR^2+vR^2) );
    aR = sqrt(gamma*pR/rR);
    HR = ( qR(4) + pR ) / rR;
    
    % First compute the Roe Averages
    RT = sqrt(rR/rL); % r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = ( HL+RT* HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-0.5*(u^2+v^2)) );
    vn = u*nx+v*ny;
    
    % Wave speed estimates
    SL=min(vnL-aL,vn-a); SLm=min(SL,0);
    SR=max(vnR+aR,vn+a); SRp=max(SR,0);
    
    % Left and Right fluxes
    FL=[rL*vnL; rL*vnL*uL+pL*nx; rL*vnL*vL+pL*ny; rL*vnL*HL];
    FR=[rR*vnR; rR*vnR*uR+pR*nx; rR*vnR*vR+pR*ny; rR*vnR*HR];
    
    % Compute strongly interacting state
    qS = (SR*qR-SL*qL+(FL-FR))/(SR-SL);
    
    % Compute the HLLE normal flux.
    nF = (SRp*FL-SLm*FR+SLm*SRp*(qR-qL))/(SRp-SLm);
        
    % Use qS and nF to compute perpedicular flux (i.e.: (q,f)->g or (q,g)->f)
    p =(nF(2)-qS(2)^2/qS(1))*nx + (nF(3)-qS(3)^2/qS(1))*ny;
    tF=[qS(3)*nx+qS(2)*ny; (qS(2)*qS(3)/qS(1))*nx+(qS(2)^2/qS(1)+p)*ny;
       (qS(3)^2/qS(1)+p)*nx+(qS(3)*qS(2)/qS(1))*ny;
       (qS(3)*(qS(4)+p/qS(1)))*nx+(qS(2)*(qS(4)+p/qS(1)))*ny];
end

function [fooN,fooS,gooE,gooW,fNo,fSo,goE,goW,sN,sS,sE,sW] = ...
    HLLE2d(C1,C2,C3,C4,gamma,dx,dy)
    % Compute HLLE2d and HLLE1d fluxes
    % [1] J. Vides, B. Nkonga, E. Audit, A simple two-dimensional extension
    % of the HLL Riemann solver for hyperbolic systems of conservation laws,
    % Journal of Computational Physics, Volume 280, 1 January 2015.
    
    % 8-states extrapolated q-values and fluxes at corner states
    qSW = C1.q;  [fSW,gSW] = Fluxes(qSW,gamma);
    qSE = C2.q;  [fSE,gSE] = Fluxes(qSE,gamma);
    qNE = C3.q;  [fNE,gNE] = Fluxes(qNE,gamma);
    qNW = C4.q;  [fNW,gNW] = Fluxes(qNW,gamma);

    % Build 2d following ideas in refernce [1]
    qSoL= C1.q + C1.dqdx*dx/2; % q_{ i ,j+1/2}^{-} from ( i , j )
    qSoR= C2.q - C2.dqdx*dx/2; % q_{ i ,j+1/2}^{+} from ( i ,j+1)
    [qSo,fSo,gSo,sSW,sSE]=HLLE1d(qSoL,qSoR,gamma,[1,0]); % mid state, flux and wave speeds

    qoWL= C1.q + C1.dqdy*dy/2; % q_{i+1/2, j }^{-} from ( i , j )
    qoWR= C3.q - C3.dqdy*dy/2; % q_{i+1/2, j }^{-} from (i+1, j )
    [qoW,goW,foW,sWS,sWN]=HLLE1d(qoWL,qoWR,gamma,[0,1]); % mid state, flux and wave speeds

    qoEL= C2.q + C2.dqdy*dy/2; % q_{ i ,j+1/2}^{-} from ( i ,j+1)
    qoER= C4.q - C4.dqdy*dy/2; % q_{ i ,j+1/2}^{-} from (i+1,j+1)
    [qoE,goE,foE,sES,sEN]=HLLE1d(qoEL,qoER,gamma,[0,1]); % mid state, flux and wave speeds

    qNoL= C3.q + C3.dqdx*dx/2; % q_{i+1/2, j }^{-} from (i+1, j )
    qNoR= C4.q - C4.dqdx*dx/2; % q_{i+1/2, j }^{-} from (i+1,j+1)
    [qNo,fNo,gNo,sNW,sNE]=HLLE1d(qNoL,qNoR,gamma,[1,0]); % mid state, flux and wave speeds

	% Verify!
    % [x,y] = meshgrid([-dx/2,0,dx/2],[-dy/2,0,dy/2]);
    % surf(x,y,zeros(3)); hold on; dt = 0.1;
    % xs = [sWS*dt,sSE*dt,sNW*dt,sNE*dt,0,0]';
    % ys = [sSW*dt,sES*dt,sWN*dt,sEN*dt,0,0]';
    % zs = [dt,dt,dt,dt,0,dt]';
    % DT = delaunayTriangulation(xs,ys,zs);
    % scatter3([sE*dt,sW*dt,0,0],[0,0,sN*dt,sS*dt],[dt,dt,dt,dt],...
    %    'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
    % tetramesh(DT); hold off; camorbit(10,0)
    
    % Verify
    % array2table([qSW,qSE,qNW,qNE],'VariableNames',{'qSW','qSE','qNW','qNE'})
    % array2table([fSW,fSE,fNW,fNE],'VariableNames',{'fSW','fSE','fNW','fNE'})
    % array2table([gSW,gSE,gNW,gNE],'VariableNames',{'gSW','gSE','gNW','gNE'})
    % array2table([qoW,qSo,qNo,qoE],'VariableNames',{'qoW','qSo','qNo','qoE'})
    % array2table([foW,fSo,fNo,foE],'VariableNames',{'foW','fSo','fNo','foE'})
    % array2table([goW,gSo,gNo,goE],'VariableNames',{'goW','gSo','gNo','goE'})        

    % Restrict certain crossings
    if (sNE < 0) && (sEN <0)        % Northeast
        if sWN > 0; sEN = 0; end
        if sSE > 0; sNE = 0; end
    end
    if (sNW < 0) && (sWN <0)        % Northwest
        if sSW > 0; sNW = 0; end
        if sEN > 0; sWN = 0; end
    end
    if (sSW < 0) && (sWS <0)        % Southwest
        if sNW > 0; sSW = 0; end
        if sES > 0; sWS = 0; end
    end
    if (sSE < 0) && (sES <0)        % Southeast
        if sNE > 0; sSE = 0; end
        if sWS > 0; sES = 0; end
    end
    
    % Fix minimun and maximun speeds
    % sN=max(sEN,sWN);	sE=max(sNE,sSE);
    % sS=max(sES,sWS);	sW=max(sNW,sSW);

    % Precompute deltas
    dq1 = sNW*sEN-sWN*sNE; df1 = sWN-sEN; dg1 = sNE-sNW;
    dq2 = sSW*sWN-sWS*sNW; df2 = sWS-sWN; dg2 = sNW-sSW;
    dq3 = sSE*sWS-sES*sSW; df3 = sES-sWS; dg3 = sSW-sSE;
    dq4 = sNE*sES-sEN*sSE; df4 = sEN-sES; dg4 = sSE-sNE;

    % Precompute c1 and c2
    c1 = dq1*dq3*(qNo-qSo) + df1*dq3*fNo - df3*dq1*fSo + dg1*dq3*gNo - dg3*dq1*gSo;
    c2 = dq4*dq2*(qoE-qoW) + df4*dq2*foE - df2*dq4*foW + dg4*dq2*goE - dg2*dq4*goW;

    % Precompute elements of inv(AII) = 1/(a*d-b*c)*[d,-b;-c,a]
    a11 = df1*dq3-df3*dq1;    a12 = dg1*dq3-dg3*dq1;
    a21 = df4*dq2-df2*dq4;    a22 = dg4*dq2-dg2*dq4;

    % Compute fluxes of the Strongly Interacting state: f** and g**
    foo=( a22*c1-a12*c2)/(a11*a22-a12*a21);
    goo=(-a21*c1+a11*c2)/(a11*a22-a12*a21);

    % Define speeds \tilde{s}_alpha for alpha \in (N,S,E,W)
    if (sES>=0) && (sWS>=0)     % Above x-axis
        sE = sSE;
        sW = sSW;
    elseif (sEN<=0) && (sWN<=0) % Below x-axis
        sE = sNE;
        sW = sNW;
    else
        sE = max(sNE,0)-max(sEN,0)*(max(sSE,0)-max(sNE,0))/(min(sES,0)-max(sEN,0));
        sW = min(sSW,0)-min(sWS,0)*(min(sNW,0)-min(sSW,0))/(max(sWN,0)-min(sWS,0));
    end
    if (sNW>=0) && (sSW>=0)     % Right of y-axis
        sN = sWN;
        sS = sWS;
    elseif (sNE<=0) && (sSE<=0) % Left of y-axis
        sN = sEN;
        sS = sES;
    else
        sN = max(sWN,0)-min(sNW,0)*(max(sWN,0)-max(sEN,0))/(min(sNW,0)-max(sNE,0));
        sS = min(sES,0)-max(sSE,0)*(min(sES,0)-min(sWS,0))/(max(sSE,0)-min(sSW,0));
    end

    % Define fluxes phiN^{~HLL2D}, phiS^{~HLL2D}, phiE^{~HLL2D} and phiW^{~HLL2D}
    sY = max(abs(sN),abs(sS));
    sX = max(abs(sE),abs(sW));
    %
    if (sW>=0) && (sS>=0)
        fooN = ((sN-sS)*foW+sS*fSW)/sN;
        fooS = fSW;
        gooE = ((sE-sW)*gSo+sW*gSW)/sE;
        gooW = gSW;
    elseif (sW>=0) && (sN<=0)
        fooN = fNW;
        fooS = ((sS-sN)*foW+sN*fNW)/sS;
        gooE = ((sE-sW)*gSo+sW*gNW)/sE;
        gooW = gNW;
    elseif (sE<=0) && (sS>=0)
        fooN = ((sN-sS)*foE+sS*fSE)/sN;
        fooS = fSE;
        gooE = gSE;
        gooW = ((sW-sE)*gSo+sE*gSE)/sW;
    elseif (sE<=0) && (sN<=0)
        fooN = fNE;
        fooS = ((sS-sN)*foE+sN*fNE)/sS;
        gooE = gNE;
        gooW = ((sW-sE)*gNo+sE*gNE)/sW;
    elseif sW>=0
        fooN = ((sY+sN)*fNW-sN*foW)/sY;
        fooS = ((sY-sS)*fSW+sS*goW)/sY;
        gooE = ((sE-sW)*goo+sW*goW)/sE;
        gooW = goW;
    elseif sE<=0
        fooN = ((sY-sN)*fNE+sN*foE)/sY;
        fooS = ((sY+sN)*fSE-sS*foE)/sY;
        gooE = goE;
        gooW = ((sW-sE)*goo+sE*goE)/sW;
    elseif sS>=0
        fooN = ((sN-sS)*foo+sS*fSo)/sN;
        fooS = fSo;
        gooE = ((sX-sE)*gSE+sE*gSo)/sX;
        gooW = ((sX+sW)*gSW-sW*gSo)/sX;
    elseif sN<=0
        fooN = fNo;
        fooS = ((sS-sN)*foo+sN*fNo)/sS;
        gooE = ((sX-sE)*gNE+sE*gNo)/sX;
        gooW = ((sX+sW)*gNW-sW*gNo)/sX;
    else
        fooN = ((sY-sN)*fNo+sN*foo)/sY;
        fooS = ((sY+sS)*fSo-sS*foo)/sY;
        gooE = ((sX-sE)*goE+sE*goo)/sX;
        gooW = ((sX+sW)*goW-sW*goo)/sX;
    end
end