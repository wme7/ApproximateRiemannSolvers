function [wave,s,amdq,apdq] = euler_roe_1D(q_l,q_r,gamma)
%     r"""
%     Roe Euler solver in 1d
%     
%     *aug_global* should contain -
%      - *gamma* - (float) Ratio of the heat capacities
%      - *gamma1* - (float) :math:`1 - \gamma`
%      - *efix* - (bool) Whether to use an entropy fix or not
%     
%     See :ref:`pyclaw_rp` for more details.
%     
%     :Version: 1.0 (2009-6-26)
%     """
    
    % Problem dimensions
    num_rp = size(q_l,2);
    num_waves = 3;
    num_eqn = 3;

    % Return values
    wave = zeros( num_eqn, num_rp, num_waves );
    s = zeros( num_waves, num_rp );
    amdq = zeros( num_eqn, num_rp );
    apdq = zeros( num_eqn, num_rp );
    
    % Solver parameters
    gamma1 = gamma-1;

    % Calculate Roe averages
    [u,a,enthalpy,~,~] = roe_averages(q_l,q_r,gamma);

    % Find eigenvector coefficients
    delta = q_r - q_l;
    a2 = gamma1 ./ a.^2 .* ((enthalpy -u.^2).*delta(1,:) + u.*delta(2,:) - delta(3,:));
    a3 = (delta(1,:) + (a-u) .* delta(1,:) - a.*a2) ./ (2*a);
    a1 = delta(1,:) - a2 - a3;
    
    % Compute the waves
    wave(1,:,1) = a1;
    wave(2,:,1) = a1 .* (u-a);
    wave(3,:,1) = a1 .* (enthalpy - u.*a);
    s(1,:) = u - a;
    
    wave(1,:,2) = a2;
    wave(2,:,2) = a2 .* u;
    wave(3,:,2) = a2 .* 0.5 .* u.^2;
    s(2,:) = u;
    
    wave(1,:,3) = a3;
    wave(2,:,3) = a3 .* (u+a);
    wave(3,:,3) = a3 .* (enthalpy + u.*a);
    s(3,:) = u + a;
    
    % Godunov update
    s_index = zeros(2,num_rp);
    for m = 1:num_eqn
        for mw = 1:num_waves
            s_index(1,:) = s(mw,:);
            amdq(m,:) = amdq(m,:) + min(s_index,[],1) .* wave(m,:,mw);
            apdq(m,:) = apdq(m,:) + max(s_index,[],1) .* wave(m,:,mw);
        end
    end
end