function res = FV_WENO5HLL1d(w,nx,dx,fluxMethod)
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

global gamma
R=3; I=R:(nx-R); % R: stencil size

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
qn  = w0n.*(2*vmm - 7*vm + 11*v)/6 ...
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
qp  = w0p.*( -umm + 5*um + 2*u  )/6 ...
	+ w1p.*( 2*um + 5*u  - up   )/6 ...
	+ w2p.*(11*u  - 7*up + 2*upp)/6;

%% Compute finite volume residual term, df/dx.
res=zeros(size(w)); flux=zeros(size(w));
for j = I % for all faces of the domain cells
    % compute flux at i+1/2
    %flux(:,j) = LFflux(qn(:,j-2),qp(:,j-2),gamma,smax);
    switch fluxMethod
        case 'ROE' % Roe
            flux(:,j) = ROEflux(qn(:,j-2),qp(:,j-2),gamma);
        case 'RUS' % Rusanov
            flux(:,j) = RUSflux(qn(:,j-2),qp(:,j-2),gamma);
        case 'HLLE' % HLLE
            flux(:,j) = HLLEflux(qn(:,j-2),qp(:,j-2),gamma);
        case 'AUSM' % AUSM
            flux(:,j) = AUSMflux(qn(:,j-2),qp(:,j-2),gamma);
        case 'HLLC' % HLLC
            flux(:,j) = HLLCflux(qn(:,j-2),qp(:,j-2),gamma);
    end
    % Flux contribution to the residual of every cell
    res(:, j ) = res(:, j ) + flux(:,j)/dx;
    res(:,j+1) = res(:,j+1) - flux(:,j)/dx;
end

% Flux contribution of the LEFT MOST FACE: left face of cell j=1.
%flux(:,3) = LFflux(qn(:,R),qp(:,R),gamma,smax);
switch fluxMethod
    case 'ROE' % Roe
        flux(:,3) = ROEflux(qn(:,3),qp(:,3),gamma);
    case 'RUS' % Rusanov
        flux(:,3) = RUSflux(qn(:,3),qp(:,3),gamma);
    case 'HLLE' % HLLE
        flux(:,3) = HLLEflux(qn(:,3),qp(:,3),gamma);
    case 'AUSM' % AUSM
        flux(:,3) = AUSMflux(qn(:,3),qp(:,3),gamma);
    case 'HLLC' % HLLC
        flux(:,3) = HLLCflux(qn(:,3),qp(:,3),gamma);
end
res(:,3) = res(:,3) - flux(:,3)/dx;
 
% Flux contribution of the RIGHT MOST FACE: right face of cell j=nx-1.
%flux(:,nx-2) = LFflux(qn(:,nx-1-2*R),qp(:,nx-1-2*R),gamma,smax);
switch fluxMethod
    case 'ROE' % Roe
        flux(:,nx-2) = ROEflux(qn(:,nx-1-2*R),qp(:,nx-1-2*R),gamma);
    case 'RUS' % Rusanov
        flux(:,nx-2) = RUSflux(qn(:,nx-1-2*R),qp(:,nx-1-2*R),gamma);
    case 'HLLE' % HLLE
        flux(:,nx-2) = HLLEflux(qn(:,nx-1-2*R),qp(:,nx-1-2*R),gamma);
    case 'AUSM' % AUSM
        flux(:,nx-2) = AUSMflux(qn(:,nx-1-2*R),qp(:,nx-1-2*R),gamma);
    case 'HLLC' % HLLC
        flux(:,nx-2) = HLLCflux(qn(:,nx-1-2*R),qp(:,nx-1-2*R),gamma);
end
res(:,nx-2) = res(:,nx-2) + flux(:,nx-2)/dx;

end % FVM WENO

function Roe = ROEflux(qL,qR,gamma)
    % Roe flux function
    %
    
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

function Rusanov = RUSflux(qL,qR,gamma)
    % Rusanov flux 
    
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

function AUSM = AUSMflux(qL,qR,gamma)
    % AUSM numerical flux
    %
    % M.-S. Liou and C. J. Steffen, A New Flux Splitting Scheme, Journal of
    % Computational Physics, 107, pp. 23-39, 1993.

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
    if (ML <= -1);
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

function HLLE = HLLEflux(qL,qR,gamma)
    % Compute HLLE flux

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
    if SLm >= 0  % Right-going supersonic flow
        HLLE = FL;
    elseif (SLm <= 0) && (SRp >= 0) % Subsonic flow
        select = 1;
        switch select
            case 1 % True HLLE function
                HLLE = ( SRp*FL - SLm*FR + SLm*SRp*(qR-qL) )/(SRp-SLm);
            case 2 % Rusanov flux ( as suggested by Toro's book )
                smax = max(abs(SLm),abs(SRp)); 
                HLLE = ( FR + FL + smax*(qL-qR) )/2; % Rusanov flux%
        end
    elseif  SRp <= 0 % Left-going supersonic flow
        HLLE = FR;
    end
end

% function HLLC = HLLCflux(qL,qR,gamma)
%     % Compute HLLC flux (Not finished!)
% 
%     % Left state
%     rL = qL(1);
%     uL = qL(2)./rL;
%     EL = qL(3)./rL;
%     pL = (gamma-1)*( qL(3) - rL*uL*uL/2 );
%     aL = sqrt(gamma*pL/rL);
%     %HL = ( qL(3) + pL )./ rL;
%     
%     % Right state
%     rR = qR(1);
%     uR = qR(2)./rR;
%     ER = qR(3)./rR;
%     pR = (gamma-1)*( qR(3) - rR*uR*uR/2 );
%     aR = sqrt(gamma*pR/rR);
%     %HR = ( qR(3) + pR )./ rR;
% 
%     % Compute guess pressure from PVRS Riemann solver
%     CUP  = 0.25*(rL + rR)*(aL + aR);
%     PPV  = max( 0 ,0.5*(pL + pR) + 0.5*(uL - uR)*CUP);
%     pmin = min(pL,pR);
%     pmax = max(pL,pR);
%     Qmax = pmax/pmin;
%     Quser= 2.0; % <--- parameter manually set! I don't like this :/
%     
%      if (Qmax <= Quser) && (pmin <= PPV) && (PPV <= pmax)
%      % Select PRVS Riemann solver
%          pM = PPV;
%          uM = 0.5*(uL + uR) + 0.5*(pL - pR)/CUP;
%       else
%          if PPV < pmin
%          % Select Two-Rarefaction Riemann solver
%             PQ  = (pL/pR)^(gamma - 1.0)/(2.0*gamma);
%             uM  = (PQ*uL/aL + uR/aR + 2/(gamma-1)*(PQ - 1.0))/(PQ/aL + 1.0/aR);
%             PTL = 1 + (gamma-1)/2.0*(uL - uM)/aL;
%             PTR = 1 + (gamma-1)/2.0*(uM - uR)/aR;
%             pM  = 0.5*(pL*PTL^(2*gamma/(gamma-1)) + pR*PTR^(2*gamma/(gamma-1)));
%          else 
%          % Use Two-Shock Riemann solver with PVRS as estimate
%             GEL = sqrt((2/(gamma+1)/rL)/((gamma-1)/(gamma+1)*pL + PPV));
%             GER = sqrt((2/(gamma+1)/rR)/((gamma-1)/(gamma+1)*pR + PPV));
%             pM  = (GEL*pL + GER*pR - (uR - uL))/(GEL + GER);
%             uM  = 0.5*(uL + uR) + 0.5*(GER*(pM - pR) - GEL*(pM - pL));
%          end
%       end
% 
%     % Estimate wave speeds: SL, SM and SR
% 	if pM<=pL; SL=uL-aL; else SL=uL-aL*sqrt(1+(gamma+1)/(2*gamma)*(pM/pL - 1)); end 
%     SM = uM;
% 	if pM<=pR; SR=uR+aR; else SR=uR+aR*sqrt(1+(gamma+1)/(2*gamma)*(pM/pR - 1)); end
% 
%     % Left and Right fluxes
%     FL=[rL.*uL; rL.*uL.^2+pL; uL.*(rL.*EL+pL)];
%     FR=[rR.*uR; rR.*uR.^2+pR; uR.*(rR.*ER+pR)];
%     
%     % Compute the HLL flux.
%     if SL >= 0  % Right-going supersonic flow
%         HLLC = FL;
%     elseif (SL <= 0) && (SR >= 0) % Subsonic flow
%         if SM > 0       % Subsonic flow to the right
%             ENEL = qL(3)/rL  + (SM - uL)*(SM + pL/(rL*(SL - uL)));
%             qsL(1,1) = rL*(SL - uL)/(SL - SM);
%             qsL(2,1) = qsL(1)*SM;
%             qsL(3,1) = qsL(1)*ENEL;
%             %HLLC = ( SR*FL - SL*FR + SL*SR*(qR-qL) )/(SR-SL);
%             HLLC = FL + SL*(qsL - qL);
%         else % SM <0    % Subsonic flow to the Left
%             ENER = qR(3)/rR + (SM - uR)*(SM + pR/(rR*(SR - uR)));
%             qsR(1,1) = rR*(SR - uR)/(SR - SM);
%             qsR(2,1) = qsR(1)*SM;
%             qsR(3,1) = qsR(1)*ENER;
%             %HLLC = ( SR*FL - SL*FR + SL*SR*(qR-qL) )/(SR-SL);
%             HLLC = FR + SR*(qsR - qR);
%         end
%     elseif  SR <= 0 % Left-going supersonic flow
%         HLLC = FR;
%     end
% end

function HLLC = HLLCflux(qL,qR,gamma)

    % Left state
    rL = qL(1);
    uL = qL(2)./rL;
    EL = qL(3)./rL;
    pL = (gamma-1)*( EL - 0.5*rL*uL*uL );
    HL = ( EL + pL )./ rL;
    aL = sqrt(gamma*pL/rL);
    
    % Right state
    rR = qR(1);
    uR = qR(2)./rR;
    ER = qR(3)./rR;
    pR = (gamma-1)*( ER - 0.5*rR*uR*uR );
    HR = ( ER + pR )./ rR;
    aR = sqrt(gamma*pR/rR);
    
    % Evaluate the two wave speeds: Einfeldt.
    RT = sqrt(rR/rL);
    u = (uL+RT*uR)/(1+RT);
    H = (HL+RT*HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-0.5*u*u) );
    
    % contravariant velocity
    Aq=Au*as+Av*bs;

  	% left and right state wave speed
	sL=min(qL-aL,Aq-Aa);
	sR=max(qR+aR,Aq+Aa);

	% middle state wave speed & pressure
    g1=dR*qR*(sR-qR)-dL*qL*(sL-qL)+pL-pR;
    g2=dR*(sR-qR)-dL*(sL-qL);
	sM=g1/g2;
    pM=dL*(qL-sL)*(qL-sM)+pL;
    
    % Interface flux
    if (sL>=0.0) % left state
        %F=FL
        HLLC(0) = Lab*(wL(0)*qL);
        HLLC(1) = Lab*(wL(1)*qL + pL*as);
        HLLC(2) = Lab*(wL(2)*qL + pL*bs);
        HLLC(3) = Lab*(wL(3)*qL + pL*qL);

        ddt = Lab*(fabs(qL)+aL);
    elseif (sL<0.0) && (0.0<=sM) % middle left
        %F=FL_star
        f1=sL-qL;
        f2=sL-sM;
        f3=f1/f2;

        wLM(0)= wL(0)*f3;
        wLM(1)= wL(1)*f3 + (pM-pL)*as/f2;
        wLM(2)= wL(2)*f3 + (pM-pL)*bs/f2;
        wLM(3)= wL(3)*f3 + (-pL*qL+pM*sM)/f2;

        HLLC(0) =  Lab*(wLM(0)*sM);
        HLLC(1) =  Lab*(wLM(1)*sM + pM*as);
        HLLC(2) =  Lab*(wLM(2)*sM + pM*bs);
        HLLC(3) =  Lab*(wLM(3)*sM + pM*sM);

        ddt = Lab*(fabs(sM)+Aa);
    elseif (sM<0.0) && (0.0<sR) % middle right
        %F=FR_star
        f1=sR-qR;
        f2=sR-sM;
        f3=f1/f2;

        wRM(0)= wR(0)*f3;
        wRM(1)= wR(1)*f3 + (pM-pR)*as/f2;
        wRM(2)= wR(2)*f3 + (pM-pR)*bs/f2;
        wRM(3)= wR(3)*f3 + (-pR*qR+pM*sM)/f2;

        HLLC(0) = Lab*(wRM(0)*sM);
        HLLC(1) = Lab*(wRM(1)*sM + pM*as);
        HLLC(2) = Lab*(wRM(2)*sM + pM*bs);
        HLLC(3) = Lab*(wRM(3)*sM + pM*sM);

        ddt = Lab*(fabs(sM)+Aa);
    elseif (sR<=0.0) % right state
        %F=FR
        HLLC(0) = Lab*(wR(0)*qR);
        HLLC(1) = Lab*(wR(1)*qR + pR*as);
        HLLC(2) = Lab*(wR(2)*qR + pR*bs);
        HLLC(3) = Lab*(wR(3)*qR + pR*qR);

        ddt = Lab*(fabs(qR)+aR);
    end
end