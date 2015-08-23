% A1
% u5568225
% u5349877

function transfer(s)
% Plots transfer function of the butterworth filter
    close all;
    c = 1;
    storer = [0; 0];
    if (s)
        for n = 10.^(1.5:0.005:1.9)
            [a, b] = afilter(n,1/(1000*n),10/n,5/n,10^(-7),45300,11300);
            loglog(n,b/a,'Marker','+','MarkerEdgeColor','blue');
            drawnow;
            hold on;
            [a2, b2] = afilter(n,1/(1000*n),10/n,5/n,10^(-7),47000,12000);
            loglog(n,b2/a2,'Marker','+','MarkerEdgeColor','red');
            drawnow;
            hold on;
        end
    else
        for n = 10.^(0:0.01:3)
            [a, b] = afilter(n,1/(1000*n),10/n,5/n,10^(-7),45300,11300);
            loglog(n,b/a,'Marker','+','MarkerEdgeColor','blue');
            drawnow;
            hold on;
            storer(1,c) = n;
            storer(2,c) = b/a;
            c = c + 1;
        end
    end
    loglog([1 1000],[1/sqrt(2) 1/sqrt(2)],'r--')
    xlabel('Input signal frequency (Hz)');
    ylabel('Output signal gain');
    legend('Output Gain', 'Cutoff value, 0.707', 'Location', 'SouthEast');



function [p_in, p_out] = afilter(omega_in,dt,t_max,trans_cutoff, C, R1, R3)
% Filter diff solution

    % Set up
    
    t = 0:dt:t_max;
    in = sin(2*pi*omega_in*t);
    l = size(t,2);
    
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
    
    [p_in,~] = findpeaks(F_in(1:M_in));
    [p_out,~] = findpeaks(F_out(1:M_out));
    
    if(isempty(p_out))
        p_out = 0;
    end
    
function [F, M, freq] = dft(u, max_freq, dt)
% DFT function
    N = size(u,2);
    F = abs(fft(u))/N;
    df = 1/(dt*N);
    freq = 0:df:max_freq;
    M = size(freq,2);
