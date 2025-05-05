%% Setting up BSA parameters.
mp.Digits(500) % Setting precision of 500 digits.
b = mp(1.08);

% transform to [1, T/dt]

% alpha = 1.2

%Experiment 1, [1, 1e+3]
% the initial precision is about 1e-14.
% M = 50;
% N = -450;
% dt = 1e-3;
% x = [];
% for i=0:2
%     x = [x,linspace(10^(i),10^(i+1),10000)];
% end
% x = x.';
% xo = [];
% for i=-3:-1
%     xo = [xo,linspace(10^(i),10^(i+1),10000)];
% end
% xo = xo.';

%Experiment 2, [1, 1e+4]
% the initial precision is about 2e-16.
% M = 50;
% N = -450;
% dt = 1e-4;
% x = [];
% for i=0:3
%     x = [x,linspace(10^(i),10^(i+1),10000)];
% end
% x = x.';
% xo = [];
% for i=-4:-1
%     xo = [xo,linspace(10^(i),10^(i+1),10000)];
% end
% xo = xo.';

%Experiment 3, [1, 1e+5]
% the initial precision is about 2e-16.
M = 70;
N = -550;
dt = 1e-5;
ddt = 1e-4;
x = [];
for i=-1:3
    x = [x,linspace(10^(i),10^(i+1),10000)];
end
x = x.';
xo = [];
for i=-5:-1
    xo = [xo,linspace(10^(i),10^(i+1),10000)];
end
xo = xo.';

%Experiment 4, [1, 1e+6]
% the initial precision is about 2e-16.
% M = 60;
% N = -640;
% dt = 1e-6;
% x = [];
% for i=0:5
%     x = [x,linspace(10^(i),10^(i+1),1000)];
% end
% x = x.';
% xo = [];
% for i=-6:-1
%     xo = [xo,linspace(10^(i),10^(i+1),10000)];
% end
% xo = xo.';
%% Preparing initial SOE.
a = 1.5;
s = b.^(N:M); % initial s
w = log(b).*(b.^((N:M)*a)).'/gamma(a); % initial w



% Testing the error of initial SOE. 
y = exp(-x*s)*w;
func = 1 ./x.^(1.5);
error_original = abs(y - func);
max(error_original)
%% 
% To emphasize the near-origin singularity, shift the SOE a little more. 
% Otherwise the accuracy near 1 can be hard to maintian.
ww = w.*exp(-0.09*s.');

% Preparing A and B. 
A = diag(-s);
B = sqrt(abs(ww));

%% Model Reduction.
n = length(s);
T = 1e4 /2;
P = B*B.';
AA = s.' + s;
EA = exp(-T*s);


% Sine all weights are positive, P=Q, which means it only needs decompose P. 
% Preparing P
% % the Classical MR
% P = P./AA;

% Time-limited Balanced Truncation; w(t)=1
% P = P.*(1 - EA.'*EA)./AA;

% Weighted Balanced Truncation
% Weight function w(t) = (r+10)^(-1/2)
K = 10;
EWA = exp(K*AA);
ggt = expint(AA*(T+K));
gg0 = expint(AA*(K));
P = P.*EWA.*(gg0-ggt);

%SVD 
[U,S,V] = svd(P);

At = U.'*A*U;
Bt = U'* B;

%% 
p = 40;
Ad = At(1:p, 1:p);
Bd = Bt(1:p);

[V,D,W] = eig(Ad);

s_nr = - diag(D);
s_nr = s_nr.'; % new s
B_mr = Bd'*V;

w_nr = B_mr.^2;
w_nr = w_nr.'; % new w
w_nr = w_nr.*exp(0.09*s_nr.'); % shift back
% filename = sprintf('%s/p_%d_mp.mat', folder, p);
% save(filename, 's_nr', 'w_nr'); % save as high precision
% 
% filename = sprintf('%s/p_%d_double.mat', folder, p);
% s_d = mp(s_nr,20);
% w_d = mp(w_nr,20);
s_d = double(s_nr);
w_d = double(w_nr);
% save(filename, "s_d","w_d");  % save as double
y1 = exp(-x*s_d)*w_d;
merror = abs(y1- 1 ./(x).^(1.5));% alpha = 1
loglog(x,merror)

%% 

% T=1;
yo = 1 ./(xo).^(1.5);

n = 50;
error = zeros(1,50);
folder = 'BSA_alpha_15_dt1e_5_T1_WBT';
if ~exist(folder, 'dir')
    mkdir(folder);
end
for p = 10:60
    Ad = At(1:p, 1:p);
    Bd = Bt(1:p);

    [V,D,W] = eig(Ad);

    s_nr = - diag(D);
    s_nr = s_nr.'; % new s
    B_mr = Bd'*V;

    w_nr = B_mr.^2;
    w_nr = w_nr.'; % new w
    w_nr = w_nr.*exp(0.9*s_nr.'); % shift back
    % filename = sprintf('%s/p_%d_mp.mat', folder, p);
    % save(filename, 's_nr', 'w_nr'); % save as high precision
    % 
    filename = sprintf('%s/p_%d_double.mat', folder, p);
    % s_d = mp(s_nr,20);
    % w_d = mp(w_nr,20);
    s_nr = s_nr/ddt;
    w_nr = w_nr/ddt^1.5;
    s_d = double(s_nr);
    w_d = double(w_nr);

    save(filename, "s_d","w_d");  % save as double
    y1 = exp(-xo*s_d)*w_d;
    merror = abs(y1- yo);
    error(p-9) = max(merror);% save maximum error
end
%% 
errorT1_wbt_k5_1 = error;
filename = sprintf('%s/error_10_60.mat', folder);
save(filename, "error");

%% 

plot(10:60, log10(errorT1_wbt_k5_1));
xlabel('p')
ylabel('log10(Maximum AbsError)')
legend('WBT, w(r) = 1/sqrt(r+5)', 'Location', 'Best');  % 'Location' 参数可以调整图例在图中的位置
title('alpha = 0.5, dt=1e-6, T = 1')
hold off;