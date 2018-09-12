function HLLC3d(wL,wR,kas,kbs,Lab,ddt,GI)
%
% void riemann_hllc_3d_(REAL wL(5), REAL wR(5), REAL n(3), REAL &Lab, REAL &ddt, REAL GI(5))
%
% wL  left state conservative variables \f$ \mathbf{W}_{L}\f$
% wR  right state conservative variables \f$ \mathbf{W}_{R}\f$
% n   normal vector \f$ \vec{n} \f$
% Lab absolute value of the differential coefficient \f$\lambda=S/\Omega\f$
% ddt spectral radius of the Jacobian matrix \f$ \rho(\mathbf{A})\f$
% GI  interface flux \f$ \mathbf{F}_{c}  \f$

% Implementation incomplete ;P

% -1>unify
absN=sqrt(n(0)*n(0)+n(1)*n(1)+n(2)*n(2));
if(fabs(absN-Lab)>1e-6)
    %   cout<<'Warning:'
    %   <<'\tabsN='<<absN
    %   <<'\tLab='<<Lab
    %   <<endl;
    %   exit(1);
end

N(0)=n(0)/absN;
N(1)=n(1)/absN;
N(2)=n(2)/absN;

% 0->basic variables for left and right states
dL = wL(0);
VxL= wL(1)/dL;
VyL= wL(2)/dL;
VzL= wL(3)/dL;
eL = wL(4);
pL = (GAMMA-1.0)*(eL-0.5*dL*(VxL*VxL+VyL*VyL+VzL*VzL));
hL = (eL+pL)/dL;
aL = sqrt(GAMMA*pL/dL);
qL = N(0)*VxL+N(1)*VyL+N(2)*VzL;

dR = wR(0);
VxR= wR(1)/dR;
VyR= wR(2)/dR;
VzR= wR(3)/dR;
eR = wR(4);
pR = (GAMMA-1.0)*(eR-0.5*dR*(VxR*VxR+VyR*VyR+VzR*VzR));
hR = (eR+pR)/dR;
aR = sqrt(GAMMA*pR/dR);
qR = N(0)*VxR+N(1)*VyR+N(2)*VzR;

% 1->Roe's average
f=sqrt(dR/dL);
AVx=(VxL+VxR*f)/(1+f);
AVy=(VyL+VyR*f)/(1+f);
AVz=(VzL+VzR*f)/(1+f);
Ah=(hL+hR*f)/(1+f);
Aa=sqrt((GAMMA-1.0)*(Ah-0.5*(AVx*AVx+AVy*AVy+AVz*AVz)));

% 2->contravariant velocity
Aq=AVx*N(0)+AVy*N(1)+AVz*N(2);

% 3->left state wave speed
sL=min(qL-aL,Aq-Aa);

% 4->right state wave speed
sR=max(qR+aR,Aq+Aa);

% 5->middle state wave speed & pressure
g1=dR*qR*(sR-qR)-dL*qL*(sL-qL)+pL-pR;
g2=dR*(sR-qR)-dL*(sL-qL);
sM=g1/g2;

pM=dL*(qL-sL)*(qL-sM)+pL;

% 6->the interface flux
% 6-1>left state
if(sL>=0.0)
    %F=FL
    
    GI(0) = Lab*(wL(0)*qL);
    GI(1) = Lab*(wL(1)*qL + pL*N(0));
    GI(2) = Lab*(wL(2)*qL + pL*N(1));
    GI(3) = Lab*(wL(3)*qL + pL*N(2));
    GI(4) = Lab*(wL(4)*qL + pL*qL);
    
    ddt=Lab*(fabs(qL)+aL);
end

% 6-2>middle left
if(sL<0.0&&0.0<=sM)
    %F=FL_star
    
    f1=sL-qL;
    f2=sL-sM;
    f3=f1/f2;
    
    wLM(0)= wL(0)*f3;
    wLM(1)= wL(1)*f3 + (pM-pL)*N(0)/f2;
    wLM(2)= wL(2)*f3 + (pM-pL)*N(1)/f2;
    wLM(3)= wL(3)*f3 + (pM-pL)*N(2)/f2;
    wLM(4)= wL(4)*f3 + (-pL*qL+pM*sM)/f2;
    
    GI(0) =  Lab*(wLM(0)*sM);
    GI(1) =  Lab*(wLM(1)*sM + pM*N(0));
    GI(2) =  Lab*(wLM(2)*sM + pM*N(1));
    GI(3) =  Lab*(wLM(3)*sM + pM*N(2));
    GI(4) =  Lab*(wLM(4)*sM + pM*sM);
    
    ddt = Lab*(fabs(sM)+Aa);
end

% 6-3>middle right
if(sM<0.0&&0.0<sR)
    %F=FR_star
    
    f1=sR-qR;
    f2=sR-sM;
    f3=f1/f2;
    
    wRM(0)= wR(0)*f3;
    wRM(1)= wR(1)*f3 + (pM-pR)*N(0)/f2;
    wRM(2)= wR(2)*f3 + (pM-pR)*N(1)/f2;
    wRM(3)= wR(3)*f3 + (pM-pR)*N(2)/f2;
    wRM(4)= wR(4)*f3 + (-pR*qR+pM*sM)/f2;
    
    GI(0) = Lab*(wRM(0)*sM        );
    GI(1) = Lab*(wRM(1)*sM + pM*N(0));
    GI(2) = Lab*(wRM(2)*sM + pM*N(1));
    GI(3) = Lab*(wRM(3)*sM + pM*N(2));
    GI(4) = Lab*(wRM(4)*sM + pM*sM);
    
    ddt = Lab*(fabs(sM)+Aa);
end

% 6-4>right state
if(sR<=0.0)
    %F=FR
    GI(0) = Lab*(wR(0)*qR);
    GI(1) = Lab*(wR(1)*qR + pR*N(0));
    GI(2) = Lab*(wR(2)*qR + pR*N(1));
    GI(3) = Lab*(wR(3)*qR + pR*N(2));
    GI(4) = Lab*(wR(4)*qR + pR*qR);
    
    ddt=Lab*(fabs(qR)+aR);
end
