%% Function prototype of the HLLE2d
% Coded by Manuel A. Diaz, NTU, 2015.11.08
clear; clc; close all;
global gamma; gamma = 1.4;

% Build a 2x2 mesh
dx=0.5; dy=0.5; [x,y]=meshgrid([0,1],[0,1]);
[r0,u0,v0,p0] = Euler_IC2d(x,y,05);
E0 = p0./((gamma-1)*r0)+0.5*(u0.^2+v0.^2);  % Total Energy
c0 = sqrt(gamma*p0./r0);                    % Speed of sound
 q = cat(3, r0, r0.*u0, r0.*v0, r0.*E0);    % initial state
 
% The corner data is obtained as
i=1; j=1;
qSW = squeeze( q( i , j ,:) );
qSE = squeeze( q( i ,j+1,:) );
qNW = squeeze( q(i+1, j ,:) );
qNE = squeeze( q(i+1,j+1,:) );

% West state
rSW = qSW(1);
uSW = qSW(2)/rSW;
vSW = qSW(3)/rSW;
pSW = (gamma-1)*( qSW(4) - rSW*(uSW^2+vSW^2)/2 );
aSW = sqrt(gamma*pSW/rSW);
HSW = ( qSW(4) + pSW ) / rSW;

% East state
rSE = qSE(1);
uSE = qSE(2)/rSE;
vSE = qSE(3)/rSE;
pSE = (gamma-1)*( qSE(4) - rSE*(uSE^2+vSE^2)/2 );
aSE = sqrt(gamma*pSE/rSE);
HSE = ( qSE(4) + pSE ) / rSE;

% South state
rNW = qNW(1);
uNW = qNW(2)/rNW;
vNW = qNW(3)/rNW;
pNW = (gamma-1)*( qNW(4) - rSW*(uNW^2+vNW^2)/2 );
aNW = sqrt(gamma*pNW/rNW);
HNW = ( qNW(4) + pNW ) / rNW;

% North state
rNE = qNE(1);
uNE = qNE(2)/rNE;
vNE = qNE(3)/rNE;
pNE = (gamma-1)*( qNE(4) - rNE*(uNE^2+vNE^2)/2 );
aNE = sqrt(gamma*pNE/rNE);
HNE = ( qNE(4) + pNE ) / rNE;




% Compute Roe Averages - SW to SE
rSroe = sqrt(rSE/rSW); 
uSroe = (uSW+rSroe*uSE)/(1+rSroe);
vSroe = (vSW+rSroe*vSE)/(1+rSroe);
HSroe = (HSW+rSroe*HSE)/(1+rSroe);
aSroe = sqrt( (gamma-1)*(HSroe-0.5*(uSroe^2+vSroe^2)) );

% Compute Roe Averages - NW to NE
rNroe = sqrt(rNE/rNW); 
uNroe = (uNW+rNroe*uNE)/(1+rNroe);
vNroe = (vNW+rNroe*vNE)/(1+rNroe);
HNroe = (HNW+rNroe*HNE)/(1+rNroe);
aNroe = sqrt( (gamma-1)*(HNroe-0.5*(uNroe^2+vNroe^2)) );

% Compute Roe Averages - SW to NW
rWroe = sqrt(rSE/rSW); 
uWroe = (uSW+rWroe*uSE)/(1+rWroe);
vWroe = (vSW+rWroe*vSE)/(1+rWroe);
HWroe = (HSW+rWroe*HSE)/(1+rWroe);
aWroe = sqrt( (gamma-1)*(HWroe-0.5*(uWroe^2+vWroe^2)) );

% Compute Roe Averages - SE to NE
rEroe = sqrt(rNE/rSE); 
uEroe = (uSE+rEroe*uNE)/(1+rEroe);
vEroe = (vSE+rEroe*vNE)/(1+rEroe);
HEroe = (HSE+rEroe*HNE)/(1+rEroe);
aEroe = sqrt( (gamma-1)*(HEroe-0.5*(uEroe^2+vEroe^2)) );




% Wave speed estimates in the S
sSW = min([ uSW-aSW, uSW+aSW, uSroe-aSroe, uSroe+aSroe ]);
sSE = max([ uSE-aSE, uSE+aSE, uSroe-aSroe, uSroe+aSroe ]);

% Wave speed estimates in the N
sNW = min([ uNW-aNW, uNW+aNW, uNroe-aNroe, uNroe+aNroe ]);
sNE = max([ uNE-aNE, uNE+aNE, uNroe-aNroe, uNroe+aNroe ]);

% Wave speed estimates in the W
sWS = min([ uSW-aSW, uSW+aSW, uWroe-aWroe , uWroe+aWroe ]);
sWN = max([ uNW-aNW, uNW+aNW, uWroe-aWroe , uWroe+aWroe ]);

% Wave speed estimates in the E
sES = min([ uSE-aSE, uSE+aSE, uEroe-aEroe, uEroe+aEroe ]);
sEN = max([ uNE-aNE, uNE+aNE, uEroe-aEroe, uEroe+aEroe ]);

% The maximum wave speed delimit the interacting region to a square domain
sS  = min(sWS,sES); 
sN  = max(sWN,sEN); 
sW  = min(sSW,sNW); 
sE  = max(sSE,sNE); 

% Velocities at cells intersections
sW_hat = sSW-sSW*(sNW-sSW)/(sWN-sWS);
sE_hat = sNE-sEN*(sSE-sNE)/(sES-sEN);
sS_hat = sES-sSE*(sES-sWS)/(sSE-sSW);
sN_hat = sWN-sNW*(sWN-sEN)/(sNW-sNE);


% Verify, Verify, Verify!
[x,y] = meshgrid([-dx/2,0,dx/2],[-dy/2,0,dy/2]);
surf(x,y,zeros(3)); hold on; dt = 0.1;
xs = [sSW*dt,sNW*dt,sNE*dt,sSE*dt,0,0]';
ys = [sWS*dt,sWN*dt,sEN*dt,sES*dt,0,0]';
zs = [dt,dt,dt,dt,0,dt]';
DT = delaunayTriangulation(xs,ys,zs);
scatter3([sNE,sNW,sSE,sSW]*dt,[sEN,sWN,sES,sWS]*dt,[dt,dt,dt,dt]);
rectangle('Position',[sW*dt sS*dt (sE-sW)*dt (sN-sS)*dt]);
scatter3([0,0,sE*dt,sW*dt],[sN*dt,sS*dt,0,0],[dt,dt,dt,dt],...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
scatter3([0,0,sE_hat*dt,sW_hat*dt],[sN_hat*dt,sS_hat*dt,0,0],[dt,dt,dt,dt],...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1]);
xlabel('x'); ylabel('y'); tetramesh(DT); hold off; view(0,90);



% Compute fluxes
fSW = [rSW*uSW; rSW*uSW*uSW + pSW; rSW*vSW*uSW; rSW*uSW*HSW];
fSE = [rSE*uSE; rSE*uSE*uSE + pSE; rSE*vSE*uSE; rSE*uSE*HSE];
fNW = [rNE*uNE; rNE*uNE*uNE + pNE; rNE*vNE*uNE; rNE*uNE*HNE];
fNE = [rNE*uNE; rNE*uNE*uNE + pNE; rNE*vNE*uNE; rNE*uNE*HNE];

gSW = [rSW*vSW; rSW*vSW*uSW; rSW*vSW*vSW + pSW; rSW*vSW*HSW];
gSE = [rSE*vSE; rSE*vSE*uSE; rSE*vSE*vSE + pSE; rSE*vSE*HSE];
gNW = [rNW*vNW; rNW*vNW*uNW; rNW*vNW*vNW + pNW; rNW*vNW*HNW];
gNE = [rNE*vNE; rNE*vNE*uNE; rNE*vNE*vNE + pNE; rNE*vNE*HNE];

% Compute the intermediate states
qSO = ( sSE*qSE - sSW*qSW + fSW-fSE )/(sSE-sSW);
qNO = ( sNE*qNE - sNW*qNW + fNW-fNE )/(sNE-sNW);
qOW = ( sWN*qNW - sWS*qSW + gSW-gNW )/(sWN-sWS);
qOE = ( sEN*qNE - sES*qSE + gSE-gNE )/(sEN-sES);

% Compute the intermediate states fluxes (normal HLLE 1d fluxes)
fSO = ( sSE*fSW - sSW*fSE + sSW*sSE*(qSE-qSW) )/(sSE-sSW);
fNO = ( sNE*fNW - sNW*fNE + sNW*sNE*(qNE-qNW) )/(sNE-sNW);
gOW = ( sWN*gSW - sWS*gNW + sWS*sWN*(qNW-qSW) )/(sWN-sWS);
gOE = ( sEN*gSE - sES*gNE + sES*sEN*(qNE-qSE) )/(sEN-sES);

% Compute the transverse intermediate fluxes (Balsara's solution)
fOW = [qOW(2);gOW(3)+(qOW(2)^2-qOW(3)^2)/qOW(1);qOW(3)*qOW(2)/qOW(1);qOW(2)*gOW(4)/qOW(3)];
fOE = [qOE(2);gOE(3)+(qOW(2)^2-qOE(3)^2)/qOE(1);qOE(3)*qOE(2)/qOE(1);qOE(2)*gOE(4)/qOE(3)];
gSO = [qSO(3);qSO(2)*qSO(3)/qSO(1);fSO(2)+(qSO(3)^2-qSO(2)^2)/qSO(1);qSO(3)*fSO(4)/qSO(2)];
gNO = [qNO(3);qNO(2)*qNO(3)/qNO(1);fNO(2)+(qNO(3)^2-qNO(2)^2)/qNO(1);qNO(3)*fNO(4)/qNO(2)];


% Area of the main quadrilateral
aOO = (dt^2/2)*((sNE-sSW)*(sWN-sES)+(sNE-sWS)*(sSE-sNW)); disp(aOO);
a22 = polyarea([sNE,sNW,sSE,sSW]*dt,[sEN,sWN,sES,sWS]*dt); disp(a22);
disp(aOO==a22)

% Strongly Interacting state q**
qOO = 1/((sNE-sSW)*(sWN-sES)+(sNE-sWS)*(sSE-sNW))*( ...
     (sWN*sNE+sSE*sEN)*qNE - (sEN*sNW+sSW*sWN)*qNW + ...
     (sES*sSW+sNW*sWN)*qSW - (sWS*sSE+sNE*sES)*qSE ...
   - sWN*fNE+sEN*fNW - sES*fSW+sWS*fSE - (sEN-sES)*fOE+(sWN-sWS)*fOW ...
   - sSE*gNE+sSW*gNW - sNW*gSW+sNE*gSE - (sNE-sNW)*gNO+(sSE-sSW)*gSO );

%% Form 0: Compute fluxes of the strongly interacting state
% Precompute deltas
dq1 = sNW*sEN-sWN*sNE; df1 = sWN-sEN; dg1 = sNE-sNW;
dq2 = sSW*sWN-sWS*sNW; df2 = sWS-sWN; dg2 = sNW-sSW;
dq3 = sSE*sWS-sES*sSW; df3 = sES-sWS; dg3 = sSW-sSE;
dq4 = sNE*sES-sEN*sSE; df4 = sEN-sES; dg4 = sSE-sNE;

% Using LSQ
b1 = dq1*(qNO-qOO) + df1*fNO + dg1*gNO;
b2 = dq2*(qOW-qOO) + df2*fOW + dg2*gOW;
b3 = dq3*(qSO-qOO) + df3*fSO + dg3*gSO;
b4 = dq4*(qOE-qOO) + df4*fOE + dg4*gOE;

% k-weights
k11 = df1*(dg2^2+dg3^2+dg4^2) - dg1*(df2*dg2+df3*dg3+df4*dg4);
k12 = df2*(dg1^2+dg3^2+dg4^2) - dg2*(df1*dg1+df3*dg3+df4*dg4);
k13 = df3*(dg1^2+dg2^2+dg4^2) - dg3*(df1*dg1+df2*dg2+df4*dg4);
k14 = df4*(dg1^2+dg2^2+dg3^2) - dg4*(df1*dg1+df2*dg2+df3*dg3);
k21 = dg1*(df2^2+df3^2+df4^2) - df1*(df2*dg2+df3*dg3+df4*dg4);
k22 = dg2*(df1^2+df3^2+df4^2) - df2*(df1*dg1+df3*dg3+df4*dg4);
k23 = dg3*(df1^2+df2^2+df4^2) - df3*(df1*dg1+df2*dg2+df4*dg4);
k24 = dg4*(df1^2+df2^2+df3^2) - df4*(df1*dg1+df2*dg2+df3*dg3);

%A = [df1,dg1;df2,dg2;df3,dg3;df4,dg4]; M=A'*A; detM=det(M);
detM = (df1*dg2-df2*dg1)^2 + (df1*dg3-df3*dg1)^2 + (df2*dg4-df4*dg2)^2 + ...
       (df3*dg2-df2*dg3)^2 + (df4*dg1-df1*dg4)^2 + (df4*dg3-df3*dg4)^2 ; % verified!

% compute fluxes of Strongly Interacting state f** and g** 
fOO = (k11*b1 + k12*b2 + k13*b3 + k14*b4)/detM;
gOO = (k21*b1 + k22*b2 + k23*b3 + k24*b4)/detM;

%% Form I
A1=[sSE+sSW-sNE-sNW,sNW+sSW-sNE-sSE;...
    sES+sWS-sEN-sWN,sWN+sWS-sEN-sES];

f00=zeros(4,1); g00=zeros(4,1);
for i=1:4
    F00=-dt^2/(4*aOO)*A1*[b1(i)-b3(i);b4(i)-b2(i)]; 
    f00(i)=F00(1); g00(i)=F00(2);
end


%% Form II
% Precompute c1 and c2
c1 = dq1*dq3*(qNO-qSO) + df1*dq3*fNO - df3*dq1*fSO + dg1*dq3*gNO - dg3*dq1*gSO;
c2 = dq4*dq2*(qOE-qOW) + df4*dq2*fOE - df2*dq4*fOW + dg4*dq2*gOE - dg2*dq4*gOW;

% Precompute elements of inv(AII) = 1/(a*d-b*c)*[d,-b;-c,a]
a11 = df1*dq3-df3*dq1;    a12 = dg1*dq3-dg3*dq1;
a21 = df4*dq2-df2*dq4;    a22 = dg4*dq2-dg2*dq4;

% Compute fluxes of the Strongly Interacting state: f** and g**
foo=( a22*c1-a12*c2)/(a11*a22-a12*a21);
goo=(-a21*c1+a11*c2)/(a11*a22-a12*a21);

% Compare solutions of methods
disp([fOO,f00,foo]);
disp([gOO,g00,goo]);
