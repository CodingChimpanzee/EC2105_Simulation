function [result] = inverse_fourier(ft, M, N)
%   INVERSE_FOURIER_TRANSFORM for total patch
%   M = w, N = h

const = 2*pi;

% for loop for m : M by M Matrix
sum1 = ones(M, M);
for u = 1 : M
    sum1(u, u) = 1;
    for m = 2 : M
        sum1(u, u) = complex(cos(const*(m-1)*(u-1)/M), sin(const*(m-1)*(u-1)/M));
    end
end


% for loop for n : N by N Matrix
sum2 = ones(N, N);
for v = 1 : N
    sum2(v, v) = 1;
    for n = 2 : N
        sum2(v, v) = complex(cos(const*(n-1)*(v-1)/N), sin(const*(n-1)*(v-1)/N));
    end
end

% Making final fourier transform: Conjugate
% Result is F itself
ft = double(ft);
result = complex(real(sum2)*real(ft), imag(sum2)*imag(ft));
result = complex(real(result)*real(sum1), imag(result)*imag(sum1));
end

