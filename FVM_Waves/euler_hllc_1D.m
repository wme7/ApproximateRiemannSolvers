function [wave,s,amdq,apdq] = euler_hllc_1D(q_l,q_r,gamma)
%     HLLC Euler solver ::
% 
%     W_1 = q_hat_l - q_l      s_1 = min(u_l-c_l,u_l+c_l,lambda_roe_1,lambda_roe_2)
%     W_2 = q_hat_r - q_hat_l  s_2 = s_m
%     W_3 = q_r - q_hat_r      s_3 = max(u_r-c_r,u_r+c_r,lambda_roe_1,lambda_roe_2)
% 
%     s_m = (p_r - p_l + rho_l*u_l*(s_l - u_l) - rho_r*u_r*(s_r - u_r))\
%           / (rho_l*(s_l-u_l) - rho_r*(s_r - u_r))
% 
%     % left middle state
%     q_hat_l(0,:) = rho_l*(s_l - u_l)/(s_l - s_m)
%     q_hat_l(1,:) = rho_l*(s_l - u_l)/(s_l - s_m)*s_m
%     q_hat_l(2,:) = rho_l*(s_l - u_l)/(s_l - s_m)\
%                 *(E_l/rho_l + (s_m - u_l)*(s_m + p_l/(rho_l*(s_l - u_l))))
%     
%     % right middle state
%     q_hat_r(0,:) = rho_r*(s_r - u_r)/(s_r - s_m)
%     q_hat_r(1,:) = rho_r*(s_r - u_r)/(s_r - s_m)*s_m
%     q_hat_r(2,:) = rho_r*(s_r - u_r)/(s_r - s_m)\
%                 *(E_r/rho_r + (s_m - u_r)*(s_m + p_r/(rho_r*(s_r - u_r))))
% 
%     *problem_data* should contain:
%     - *gamma* - (float) Ratio of specific heat capacities
%     - *gamma1* - (float) :math:`\gamma - 1`
% 
%     :Version 1.0 (2015-11-18)

    % Problem dimensions
    num_rp = size(q_l,2);
    num_waves = 3;
    num_eqn = 3;

    % Return values
    wave = zeros( num_eqn, num_rp, num_waves);
    s = zeros( num_waves, num_rp );
    amdq = zeros( num_eqn, num_rp );
    apdq = zeros( num_eqn, num_rp );
    
    % Calculate Roe averages, right and left speeds
    [u,a,~,p_l,p_r] = roe_averages(q_l,q_r,gamma);
    rho_r = q_r(1,:);
    rho_l = q_l(1,:);
    E_r = q_r(3,:);
    E_l = q_l(3,:);
    H_r = (E_r + p_r)./rho_r;
    H_l = (E_l + p_l)./rho_l;
    u_r = q_r(2,:)./rho_r;
    u_l = q_l(2,:)./rho_l;
    a_r = sqrt((gamma-1)*(H_r - 0.5*u_r.^2));
    a_l = sqrt((gamma-1)*(H_l - 0.5*u_l.^2));

    % Compute Einfeldt speeds
    s_index = zeros(4,num_rp);
    s_index(1,:) = u + a;
    s_index(2,:) = u - a;
    s_index(3,:) = u_l + a_l;
    s_index(4,:) = u_l - a_l;
    s(1,:) = min(s_index,[],1);
    s_index(3,:) = u_r + a_r;
    s_index(4,:) = u_r - a_r;
    s(3,:) = max(s_index,[],1);

    % left and right speeds
    s_l = s(1,:);
    s_r = s(3,:);
    
    % middle speed
    s_m = (p_r-p_l + rho_l.*u_l.*(s_l-u_l) - rho_r.*u_r.*(s_r-u_r)) ...
          ./ (rho_l.*(s_l-u_l) - rho_r.*(s_r-u_r));
    s(2,:) = s_m;
    
    % left middle states
    q_hat_l = zeros(num_eqn,num_rp);
    q_hat_l(1,:) = rho_l.*(s_l - u_l)./(s_l - s_m);
    q_hat_l(2,:) = rho_l.*(s_l - u_l)./(s_l - s_m).*s_m;
    q_hat_l(3,:) = rho_l.*(s_l - u_l)./(s_l - s_m) ...
                .*(E_l./rho_l + (s_m - u_l).*(s_m + p_l./(rho_l.*(s_l - u_l))));
    
    % right middle state
    q_hat_r = zeros(num_eqn,num_rp);
    q_hat_r(1,:) = rho_r.*(s_r - u_r)./(s_r - s_m);
    q_hat_r(2,:) = rho_r.*(s_r - u_r)./(s_r - s_m).*s_m;
    q_hat_r(3,:) = rho_r.*(s_r - u_r)./(s_r - s_m) ...
                .*(E_r./rho_r + (s_m - u_r).*(s_m + p_r./(rho_r.*(s_r - u_r))));

    % Compute each family of waves
    wave(:,:,1) = q_hat_l - q_l;
    wave(:,:,2) = q_hat_r - q_hat_l;
    wave(:,:,3) = q_r - q_hat_r;
    
    % Compute variations
    s_index = zeros(2,num_rp);
    for m = 1:num_eqn
        for mw = 1:num_waves
            s_index(1,:) = s(mw,:);
            amdq(m,:) = amdq(m,:) + min(s_index,[],1) .* wave(m,:,mw);
            apdq(m,:) = apdq(m,:) + max(s_index,[],1) .* wave(m,:,mw);
        end
    end
    
end % Euler_HLLC_1d
