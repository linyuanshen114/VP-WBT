%% 
N = 1000;
x = linspace(0,10,N);
% M = -700;
% n=2000;
% b = 1.01;
% num = M:0;
L = 1;
% Compute n-point Gauss-Legendre nodes and weights on [a,b]
% [x0, w0] = lgwt(n, 0, 10);
% I = sum(w .* f(x));

w = 2*L/sqrt(pi)*10/N*ones(N,1);
s = L^2*x.^2;
s = s;
w = w.';

x = x.';
xx =x.^2;

sog = exp(-xx*s)*w.';
error = abs(sog-erf(L*(x+1e-16))./(x+1e-16));

plot(x,log10(error))