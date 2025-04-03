%% Setting up BSA parameters.
mp.Digits(400) % Setting precision of 400 digits.
b = mp(1.1);

% alpha = 1
% the initial precision is about 2e-17.
M = 40;
N = -436;

% alpha = 0.5
% the initial precision is about 5e-17.
% M = 40;
% N = -859;

% alpha = 0.25
% the initial precision is about 5e-17.
% b = mp(1.3);
% M = 15;
% N = -600;

% alpha = 2
% the initial precision is about 8e-17.
% M = 40;
% N = -659;

%% Preparing initial SOE.
% alpha = 1
s = b.^(N:M); % initial s
w = log(b).*s.'; % initial w

% alpha = 0.5
% s = b.^(N:M);
% w = log(b).*(b.^((N:M)/2)).'/gamma(0.5);

% alpha = 0.25
% s = b.^(N:M);
% w = log(b).*(b.^((N:M)/4)).'/gamma(0.25);

% alpha = 2
% s = b.^(N:M);
% w = log(b).*(b.^((N:M)*2)).'/gamma(2);

x = mp(1:0.01:1024);
x=x.';

% Testing the error of initial SOE. 
%y = exp(-x*s)*w;
%max(abs(y - 1 ./x))

% To emphasize the near-origin singularity, shift the SOE a little more. 
% Otherwise the accuracy near 1 can be hard to maintian.
ww = w.*exp(-0.9*s.');

% Preparing A and B. 
A = diag(-s);
B = sqrt(abs(ww));

%% Model Reduction.
n = length(s);
T = 1024/2;
P = B*B.';
AA = s.' + s;
EA = exp(-T*s);
EWA = exp(10*AA);

% Sine all weights are positive, P=Q, which means it only needs decompose P. 
% Preparing P
% the Classical MR
% P = P./AA;

% Time-limited Balanced Truncation; w(t)=1
% P = P.*(1 - EA.'*EA)./AA;

% Weighted Balanced Truncation
% Weight function w(t) = (r+10)^(-1/2)
ggt = expint(AA*(T+10));
gg0 = expint(AA*(10));
P = P.*EWA.*(gg0-ggt);

%SVD 
[U,S,V] = svd(P);

At = U.'*A*U;
Bt = U'* B;

%% Consturcting reduced SOE

n = 12;
error = zeros(1,12);
folder = 'BSA_alpha_1_WBT';
if ~exist(folder, 'dir')
    mkdir(folder);
end
for p = 24:35
    Ad = At(1:p, 1:p);
    Bd = Bt(1:p);

    [V,D,W] = eig(Ad);

    s_nr = - diag(D);
    s_nr = s_nr.'; % new s
    B_mr = Bd'*V;

    w_nr = B_mr.^2;
    w_nr = w_nr.'; % new w
    w_nr = w_nr.*exp(0.9*s_nr.'); % shift back
    filename = sprintf('%s/p_%d_mp.mat', folder, p);
    save(filename, 's_nr', 'w_nr'); % save as high precision

    filename = sprintf('%s/p_%d_double.mat', folder, p);
    s_d = double(s_nr);
    w_d = double(w_nr);
    save(filename, "s_d","w_d");  % save as double

    x = 1:0.01:1024;
    x = x.';
    y1 = exp(-x*s_nr)*w_nr;
    merror = abs(y1- 1 ./(x));% alpha = 1
    error(p-23) = max(merror);% save maximum error
end

