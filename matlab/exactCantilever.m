function f = exactCantilever(L,E,rho,A,I,mtip,n)
%EXACTCANTILEVER  Exact natural frequencies (Hz) of a uniform Euler-Bernoulli
%   cantilever with a point mass mtip at the free end (first n modes).
%   Frequency equation, with a = mtip/(rho*A*L) and b = beta*L:
%   1 + cos(b)cosh(b) + a*b*(cos(b)sinh(b) - sin(b)cosh(b)) = 0
a = mtip/(rho*A*L);
F = @(b) 1 + cos(b).*cosh(b) + a*b.*(cos(b).*sinh(b) - sin(b).*cosh(b));
x = linspace(0.05, 16, 20000);  v = F(x);  b = [];
for i = 1:numel(x)-1
    if v(i)*v(i+1) < 0
        b(end+1) = fzero(F, [x(i) x(i+1)], optimset('TolX',1e-15)); %#ok<AGROW>
    end
    if numel(b) == n, break; end
end
f = b.^2 * sqrt(E*I/(rho*A*L^4)) / (2*pi);
end
