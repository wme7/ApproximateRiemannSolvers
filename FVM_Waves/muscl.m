function [qL,qR] = muscl(q,limiter)

N = size(q,2);
dq=zeros(3,N); 
qL=zeros(3,N); 
qR=zeros(3,N);

% Compute and limit slopes
    for i = 1:3
        for j = 2:N-1 % for all internal faces
            switch limiter
                case 'MC' % MC limiter
                    % Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}
                    dqR = 2*(q(i,j+1) - q(i,j));
                    dqL = 2*(q(i,j) - q(i,j-1));
                    dqC = (q(i,j+1) - q(i,j-1))/(2);
                    dq(i,j) = minmod([dqR,dqL,dqC]);
                case 'MM' % Minmod limiter
                    % Find dq_j = minmod{fwd diff, bwd diff}
                    dqR = (q(i,j+1) - q(i,j));
                    dqL = (q(i,j) - q(i,j-1));
                    dq(i,j) = minmod([dqR,dqL]);
                case 'VA' % Van Albada limiter
                    dqR = (q(i,j+1) - q(i,j));
                    dqL = (q(i,j) - q(i,j-1));
                    dq(i,j) = vanalbada(dqR,dqL,dx);
            end
        end
    end

    % Left and Right extrapolated q-values at the boundary j+1/2
    for j = 2:N-2 % for the domain cells
        qL(:,j) = q(:, j ) + dq(:, j )/2;	% q_{j+1/2}^{+} from j
        qR(:,j) = q(:,j+1) - dq(:,j+1)/2;	% q_{j+1/2}^{-} from j+1
    end
    
end

function mm = minmod(v)
    % Using Harten's generalized definition
    % minmod: zero if opposite sign, otherwise the one of smaller magnitude.
    %m=size(v,1); mm=zeros(size(v,2),1); s=sum(sign(v),2)/m; ids=find(abs(s)==1);
    %if(~isempty(ids)); mm(ids)=s(ids).*min(abs(v(ids,:)),[],2); end
    s = sum(sign(v))/numel(v); 
    if abs(s)==1; mm = s*min(abs(v(:))); else, mm=0; end
end

function va = vanalbada(da,db,h)
    % Van Albada Slope Limiter Function
    % vanAlbada: extend the simetric formulation of the van leer limiter
    eps2=(0.3*h)^3; 
    va=0.5*(sign(da)*sign(db)+1)*((db^2+eps2)*da+(da^2+eps2)*db)/(da^2+db^2+2*eps2);
end