close all;
clear all;
p = 0.001;
freqs = 49:p:51;
m = size(freqs);
x = zeros(m);
counter = 1;
for n = 49:p:51
    val = log10(n);
    [pin, pout] = pfilter(10^val, 1/(10000*val), 10/val, 5/val);
    x(counter) = max(pout)/max(pin);
    counter = counter + 1;
end
index = find(x==min(x));
freqs(70)
