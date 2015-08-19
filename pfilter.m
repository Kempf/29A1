function [p_in, p_out] = pfilter(omega_in,dt,t_max,trans_cutoff)
% Filter diff solution

    % Set up
    
    t = 0:dt:t_max;
    in = cos(2*pi*omega_in*t);
    l = size(t,2);
    
    C = 10^(-7);
    R1 = 45300;
    R3 = 11300;
    
    A = [-1/(C*R3) 1/(C*R1) -1/(C*R1); 0 -2/(C*R1) 1/(C*R1); -1/(C*R3) 2/(R1*C) -2/(R1*C)];
    X = zeros(3,l);
    B = [0 1; 1/(C*R1) 0; 0 1];
    
    % Euler's
    
    for n = 1:l-1
        X(:,n+1) = X(:,n)+dt*(A*X(:,n)+B*[in(n); (in(n+1)-in(n))/dt]);
    end
    
    out = X(3,:);
    
    % DFT
    cut = round(l/t_max*trans_cutoff);
    if (cut == 0)
        cut = 1;
    end
    
    [F_out, M_out, ~] = dft(out(cut:end),omega_in*2, dt);
    [F_in, M_in, ~] = dft(in,omega_in*2, dt);
    
    [p_in,~] = findpeaks(F_in(1:M_in), 'threshold', 0);
    [p_out,~] = findpeaks(F_out(1:M_out), 'threshold', 0);
    
    if(isempty(p_out))
        p_out = 0;
    end
    
