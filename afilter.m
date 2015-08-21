% A1
% u5568225
% u5349877

% Filter diff solution
function out = afilter(omega_in,dt,t_max,trans_cutoff)
    
    % Set up
    
    t = 0:dt:t_max;
    in = sin(2*pi*omega_in*t);
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
    
    cut = round(trans_cutoff/t_max*l);
    if (cut == 0)
        cut = 1;
    end
    u = out(cut:end);
    N = size(u,2);
    %max_freq = 1/dt/2;
    max_freq = omega_in*2;
    % FFT
    F = abs(fft(u))/N;
    % frequency scale
    df = 1/(dt*N);
    freq = 0:df:max_freq;
    M = size(freq,2);
        
    % Plot
    close all;
    hold off;
    plot(t,in,'b');
    hold on;
    plot(t,out,'r');
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    title('Input vs. Output');
    hold off;
    
    figure;
    plot(freq,F(1:M));
    set(gca, 'XTick', 0:(max_freq/10):max_freq);
    title('DFT');
    xlabel('Frequency (Hz)');
    ylabel('|Y(f)|');
    
