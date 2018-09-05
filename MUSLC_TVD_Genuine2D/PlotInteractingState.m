%% Plot the Strongly Interacting state q**
% by Manuel Diaz

clear; clc; close all;

load('q_corners.mat'); dx=0.5; dy=0.5;

% % Build 2d following ideas in refernce [1]
% qSoL= cell( i , j ).q + cell( i , j ).dqdx*dx/2; % q_{ i ,j+1/2}^{-} from ( i , j )
% qSoR= cell( i ,j+1).q - cell( i ,j+1).dqdx*dx/2; % q_{ i ,j+1/2}^{+} from ( i ,j+1)
% [qSo,fSo,gSo,sSW,sSE]=HLLE1d(qSoL,qSoR,gamma,[1,0]); % mid state, flux and wave speeds
% 
% qoWL= cell( i , j ).q + cell( i,  j ).dqdy*dy/2; % q_{i+1/2, j }^{-} from ( i , j )
% qoWR= cell(i+1, j ).q - cell(i+1, j ).dqdy*dy/2; % q_{i+1/2, j }^{-} from (i+1, j )
% [qoW,goW,foW,sWS,sWN]=HLLE1d(qoWL,qoWR,gamma,[0,1]); % mid state, flux and wave speeds
% 
% qoEL= cell( i ,j+1).q + cell( i ,j+1).dqdy*dy/2; % q_{ i ,j+1/2}^{-} from ( i ,j+1)
% qoER= cell(i+1,j+1).q - cell(i+1,j+1).dqdy*dy/2; % q_{ i ,j+1/2}^{-} from (i+1,j+1)
% [qoE,goE,foE,sES,sEN]=HLLE1d(qoEL,qoER,gamma,[0,1]); % mid state, flux and wave speeds
% 
% qNoL= cell(i+1, j ).q + cell(i+1, j ).dqdx*dx/2; % q_{i+1/2, j }^{-} from (i+1, j )
% qNoR= cell(i+1,j+1).q - cell(i+1,j+1).dqdx*dx/2; % q_{i+1/2, j }^{-} from (i+1,j+1)
% [qNo,fNo,gNo,sNW,sNE]=HLLE1d(qNoL,qNoR,gamma,[1,0]); % mid state, flux and wave speeds
% 
% % verify
% array2table([qSW,qSE,qNW,qNE],'VariableNames',{'qSW','qSE','qNW','qNE'})
% array2table([fSW,fSE,fNW,fNE],'VariableNames',{'fSW','fSE','fNW','fNE'})
% array2table([gSW,gSE,gNW,gNE],'VariableNames',{'gSW','gSE','gNW','gNE'})
% 
% array2table([qoW,qSo,qNo,qoE],'VariableNames',{'qoW','qSo','qNo','qoE'})
% array2table([foW,fSo,fNo,foE],'VariableNames',{'foW','fSo','fNo','foE'})
% array2table([goW,gSo,gNo,goE],'VariableNames',{'goW','gSo','gNo','goE'})

% Define speeds \tilde{s}_\alpha for alpha \in (N,S,E,W)
if (sES>=0) && (sWS>=0)     % Speeds Above x-axis
    sE = sSE;
    sW = sSW;
elseif (sEN<=0) && (sWN<=0) % Speeds Below x-axis
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

% Verify!
[x,y] = meshgrid([-dx/2,0,dx/2],[-dy/2,0,dy/2]);
surf(x,y,zeros(3)); hold on; dt = 0.1;
xs = [sWS*dt,sSE*dt,sNW*dt,sNE*dt,0,0]';
ys = [sSW*dt,sES*dt,sWN*dt,sEN*dt,0,0]';
zs = [dt,dt,dt,dt,0,dt]';
DT = delaunayTriangulation(xs,ys,zs);
scatter3([sNE,sNW,sSE,sSW]*dt,[sEN,sWN,sES,sWS]*dt,[dt,dt,dt,dt]);
scatter3([sE*dt,sW*dt,0,0],[0,0,sN*dt,sS*dt],[dt,dt,dt,dt],...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
tetramesh(DT); hold off; camorbit(10,0);
    
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



% Pseudocode 3 in ref.[1], it is meant to design phix^{~HLL2D} &
% phiy^{~HLL2D} according with the author
%     if (sW>=0) && (sS>=0)
%         phi2x = ((sN-sS)*f_W + sS*fSW)/sN;
%         phi2y = ((sE-sW)*gS_ + sW*gSW)/sE;
%     elseif (sW>=0) && (sN<=0)
%         phi2x = ((sS-sN)*f_W + sN*fNW)/sS;
%         phi2y = ((sE-sW)*gS_ + sW*gNW)/sE;
%     elseif (sE<=0) && (sS>=0)
%         phi2x = ((sN-sS)*f_E + sS*fSE)/sN;
%         phi2y = ((sW-sE)*gS_ + sE*gSE)/sW;
%     elseif (sE<=0) && (sN<=0)
%         phi2x = ((sS-sN)*f_E + sN*fNE)/sS;
%         phi2y = ((sW-sE)*gN_ + sE*gNE)/sW;
%     elseif sW>=0
%         phi2x = f_W;
%         phi2y = ((sE-sW)*g__ + sW*g_W)/sE;
%     elseif sE<=0
%         phi2x = f_E;
%         phi2y = ((sW-sE)*g__ + sE*g_E)/sW;
%     elseif sS>=0
%         phi2x = ((sN-sS)*f__ + sS*fS_)/sN;
%         phi2y = gS_;
%     elseif sN<=0
%         phi2x = ((sS-sN)*f__ + sN*fN_)/sS;
%         phi2y = gN_;
%     else
%         phi2x = f__;
%         phi2y = g__;
%     end

%% DO NOT ERASE!

% % Using mid state and normal flux compute
%     if prod(abs(normal)==[1,0]);	% given f: x-flux compute g: y-flux 
%         p=nF(2)-qS(2)^2/qS(1); 
%         tF=[qS(3); qS(2)*qS(3)/qS(1); qS(3)^2/qS(1)+p; qS(3)*(qS(4)+p/qS(1))];
%         %tF = [qS(3); qS(2)*qS(3)/qS(1); nF(2)+(qS(3)^2-qS(2)^2)/qS(1); qS(3)*qS(4)/qS(2)];
%     else                            % given g: y-flux compute f: x-flux
%         p=nF(3)-qS(3)^2/qS(1); 
%         tF=[qS(2); qS(2)^2/qS(1)+p; qS(3)*qS(2)/qS(1); qS(2)*(qS(4)+p/qS(1))];
%         %tF = [qS(2); nF(3)+(qS(2)^2-qS(3)^2)/qS(1); qS(3)*qS(2)/qS(1); qS(2)*nF(4)/qS(3)];
%     end
