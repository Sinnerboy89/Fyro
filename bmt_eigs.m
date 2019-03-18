function [V, D, b] = bmt_eigs(fir)
%BMT_EIGS

s = coeffs(fir);
b = s.Numerator;
n = length(b) - 1;

% Hankel matrix
Hank = zeros(n, n);
for i=1:n
    Hank(:, i) = [b(i+1:end)'; zeros(i-1, 1)];
end

% eigens of Hankel matrix
[V, D] = eigs(Hank, n);

end

