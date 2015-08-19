% CLAB2
% u5568225

function [F, M, freq] = dft(u, max_freq, dt)
% DFT function
    N = size(u,2);
    F = abs(fft(u))/N;
    df = 1/(dt*N);
    freq = 0:df:max_freq;
    M = size(freq,2);
