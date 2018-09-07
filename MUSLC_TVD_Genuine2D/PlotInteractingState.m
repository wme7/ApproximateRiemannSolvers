%% Plot the Strongly Interacting state q**
% by Manuel Diaz
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
qSW = q( i , j ,:);
qSE = q( i ,j+1,:);
qNW = q(i+1, j ,:);
qNE = q(i+1,j+1,:);

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
xlabel('x'); ylabel('y'); tetramesh(DT); hold off; view(0,90);



% Area of the main quadrilateral
aoo = (dt^2/2)*((sNE-sSW)*(sWN-sES)+(sNE-sWS)*(sSE-sNW)); disp(aoo);
a22 = polyarea([sNE,sNW,sSE,sSW]*dt,[sEN,sWN,sES,sWS]*dt); disp(a22);
disp(aoo==a22)


% Strongly Interacting state q**
qoo = 1/((sNE-sSW)*(sWN-sES)+(sNE-sWS)*(sSE-sNW))*( ...
     (sWN*sNE+sSE*sEN)*qNE - (sEN*sNW+sSW*sWN)*qNW + ...
     (sES*sSW+sNW*sWN)*qSW - (sWS*sSE+sNE*sES)*qSE ...
   - sWN*fNE+sEN*fNW - sES*fSW+sWS*fSE - (sEN-sES)*foE+(sWN-sWS)*foW ...
   - sSE*gNE+sSW*gNW - sNW*gSW+sNE*gSE - (sNE-sNW)*gNo+(sSE-sSW)*gSo );



% Precompute deltas
dq1 = sNW*sEN-sWN*sNE; df1 = sWN-sEN; dg1 = sNE-sNW;
dq2 = sSW*sWN-sWS*sNW; df2 = sWS-sWN; dg2 = sNW-sSW;
dq3 = sSE*sWS-sES*sSW; df3 = sES-sWS; dg3 = sSW-sSE;
dq4 = sNE*sES-sEN*sSE; df4 = sEN-sES; dg4 = sSE-sNE;

%% USING LSQ
b1 = dq1*(qNo-qoo) + df1*fNo + dg1*gNo;
b2 = dq2*(qoW-qoo) + df2*foW + dg2*goW;
b3 = dq3*(qSo-qoo) + df3*fSo + dg3*gSo;
b4 = dq4*(qoE-qoo) + df4*foE + dg4*goE;

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
f00 = (k11*b1 + k12*b2 + k13*b3 + k14*b4)/detM;
g00 = (k21*b1 + k22*b2 + k23*b3 + k24*b4)/detM;

%% Form I



%% Form II
% Precompute c1 and c2
c1 = dq1*dq3*(qNo-qSo) + df1*dq3*fNo - df3*dq1*fSo + dg1*dq3*gNo - dg3*dq1*gSo;
c2 = dq4*dq2*(qoE-qoW) + df4*dq2*foE - df2*dq4*foW + dg4*dq2*goE - dg2*dq4*goW;

% Precompute elements of inv(AII) = 1/(a*d-b*c)*[d,-b;-c,a]
a11 = df1*dq3-df3*dq1;    a12 = dg1*dq3-dg3*dq1;
a21 = df4*dq2-df2*dq4;    a22 = dg4*dq2-dg2*dq4;

% Compute fluxes of the Strongly Interacting state: f** and g**
foo=( a22*c1-a12*c2)/(a11*a22-a12*a21);
goo=(-a21*c1+a11*c2)/(a11*a22-a12*a21);
