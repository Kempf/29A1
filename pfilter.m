function [ p_in, p_out ] = pfilter ( omega_in,dt,t_max,trans_cutoff )
% Filter ODE Euler's solution

    % Input and time scale
    
    t = 0:dt:t_max;
    in = cos(2*pi*omega_in*t);
    l = size(t,2);
    
    % ODE set up
    
    C = 10^(-7);
    R1 = 45300;
    R3 = 11300;
    
    A = [-1/(C*R3) 1/(C*R1) -1/(C*R1); ...
    0 -2/(C*R1) 1/(C*R1); -1/(C*R3) 2/(R1*C) -2/(R1*C)];
    X = zeros(3,l);
    B = [0 1; 1/(C*R1) 0; 0 1];
    
    % Euler's solution
    
    for n = 1:l-1
        X(:,n+1) = X(:,n)+dt*(A*X(:,n)+B*[in(n); (in(n+1)-in(n))/dt]);
    end
    
    out = X(3,:);
    
    % Transient cutoff
    
    cut = round(l/t_max*trans_cutoff);
    if (cut == 0)
        cut = 1;
    end
    
    % DFT
    
    function [F, M, freq] = dft(u, max_freq, dt)
        % DFT helper function
        N = size(u,2);
        F = abs(fft(u))/N;
        df = 1/(dt*N);
        freq = 0:df:max_freq;
        M = size(freq,2);
    end
    
    [F_out, M_out, ~] = dft(out(cut:end),omega_in*2, dt);
    [F_in, M_in, ~] = dft(in,omega_in*2, dt);
    
    % Finding peaks
    
    [p_in,~] = findpeaks(F_in(1:M_in));
    [p_out,~] = findpeaks(F_out(1:M_out));
    
    % Set gain to 0 if no peak was found
    
    if(isempty(p_out))
        p_out = 0;
    end
end