% load("F:\VP-OMR\src\HSVD\s_32.mat")
% load("F:\VP-OMR\src\HSVD\w_32.mat")
clc
clear
format long
global d
T = 10;         % 计算区间 [0, T]
N = 7500;         % 采样参数 (总采样点为 2N+1)
% tol = 1e-10;    % 奇异值截断阈值
r = 16;
L = 100;
shift = 0;
% d = 100;
% mp(d);

fHandle = @(t) erf(L*(t+1e-16))./(t+1e-16);  
fTest = @(t) erf(L*(t+1e-16))./(t+1e-16);
%% 

% n1 = 1; n2 = 30;
% maxerror_list = n1:n2;
% tGrid =linspace(0, T,10000);    % 在 [0,T] 上取 1000 个点进行评估
% fExact = fTest(tGrid);  
% for r= n1:n2
%     [alpha, c, r] = HSVD(fTest,shift, N, T, r);
%          % 真值
%     alpha = double(alpha);
%     c = double(c);
%     % 计算近似值: fApprox(t) = sum_{k=1}^{r} c_k * exp(alpha_k*t)
%     fApprox = mp(zeros(size(tGrid)),d);
%     for k = 1:r
%         fApprox = fApprox + c(k) * exp(alpha(k) * tGrid);
%     end
% 
%     absError = abs(fExact - real(fApprox));
%     relError = abs(fExact - real(fApprox)) ./ abs(fExact);
% 
%     % 计算在整个区间上的最大相对误差
%     maxRelError = max(relError);
%     maxAbsError = max(absError);
%     maxerror_list(r-n1+1) = maxAbsError;
% end
% plot(n1:n2,log10(maxerror_list));
% title(sprintf('Lambda = %.3f, shift = %s, Sample num = %.3f',L , string(shift), 2*N+1))
% xlabel('p')
% ylabel('log10(Maximum AbsError)')

%%
n1 = 1; n2 = 36;
maxerror_list = n1:n2;
tGrid =linspace(0, T,100000);    % 在 [0,T] 上取 1000 个点进行评估
fExact = fTest(tGrid);  
tSamples = linspace(-shift, T, 2*N + 1).';
tSamples_pos = tSamples(tSamples >= 0);
% ySamples = fHandle(tSamples);
ySamples_pos = fHandle(tSamples_pos);
H = tSamples + tSamples.';
H = H(1:N+1,1:N+1);
H = fHandle(H);

% 3. 对 Hankel 矩阵做 SVD, 并根据 tol 截断
% --------------------------------------
[Uh, Sh, V] = svd(H, 'econ');      % H = U*S*V',  'econ'可节省空间
sVals = diag(Sh);                 % 奇异值向量
sMax = sVals(1);                 % 最大奇异值
    % idx = find(sVals > tol * sMax);  % 保留的奇异值下标
    % r = length(idx);
for r= n1:n2
    U_r = Uh(:, 1:r);
    
    if r == 0
        % 如果所有奇异值都很小, 说明函数非常接近 0
        alpha = [];
        c = [];
        return;
    end
    
    % 4. 构造 U^(+) 和 U^(-)，并解特征值问题
    % --------------------------------------
    % 令 U^(+)=U_r(1:N,:), U^(-)=U_r(2:N+1,:)
    Uplus = U_r(1:N, :);      % 大小: N x r
    Uminus = U_r(2:N+1, :);   % 大小: N x r
    
    % 矩阵 pencil 方法: Z = pinv(Uplus) * Uminus
    Z = pinv(Uplus) * Uminus; % 尺寸: r x r
    
    % 解特征值问题
    [EigVec, EigValMat] = eig(Z);  % Z*w = lambda*w
    lambda = diag(EigValMat);     
    
    dt = (T + shift)/ (2*N);
    alpha = (1/dt) * log(lambda);
    Vmat = exp(tSamples_pos*(alpha.'));
    
    % 求解 c, 使得 Vmat * c ~ ySamples
    c = pinv(Vmat)*ySamples_pos; 
    alpha = double(alpha);
    c = double(c);
    % 计算近似值: fApprox(t) = sum_{k=1}^{r} c_k * exp(alpha_k*t)
    fApprox = zeros(size(tGrid));
    for k = 1:r
        fApprox = fApprox + c(k) * exp(alpha(k) * tGrid);
    end

    absError = abs(fExact - real(fApprox));
    relError = abs(fExact - real(fApprox)) ./ abs(fExact);

    % 计算在整个区间上的最大相对误差
    maxRelError = max(relError);
    maxAbsError = max(absError);
    maxerror_list(r-n1+1) = maxAbsError;
end
%% 

plot(n1:n2, log10(maxerror_list))

%% 
% r=n2;
% [alpha, c, r] = HSVD(fTest,shift, N, T, r); 
r = 27;
U_r = Uh(:, 1:r);

if r == 0
    % 如果所有奇异值都很小, 说明函数非常接近 0
    alpha = [];
    c = [];
    return;
end

% 4. 构造 U^(+) 和 U^(-)，并解特征值问题
% --------------------------------------
% 令 U^(+)=U_r(1:N,:), U^(-)=U_r(2:N+1,:)
Uplus = U_r(1:N, :);      % 大小: N x r
Uminus = U_r(2:N+1, :);   % 大小: N x r

% 矩阵 pencil 方法: Z = pinv(Uplus) * Uminus
Z = pinv(Uplus) * Uminus; % 尺寸: r x r

% 解特征值问题
[EigVec, EigValMat] = eig(Z);  % Z*w = lambda*w
lambda = diag(EigValMat);     

dt = T / (2*N);
alpha = (1/dt) * log(lambda);
Vmat = exp(tSamples_pos*(alpha.'));

% 求解 c, 使得 Vmat * c ~ ySamples
c = pinv(Vmat)*ySamples_pos; 
alpha = double(alpha);
c = double(c);
fApprox = zeros(size(tGrid));
for k = 1:r
    fApprox = fApprox + c(k) * exp(alpha(k) * tGrid);
end

absError = abs(fExact - real(fApprox));
% relError = abs(fExact - real(fApprox)) ./ abs(fExact);
plot(tGrid,log10(absError))

%% 
r = 30;
U_r = Uh(:, 1:r);

if r == 0
    % 如果所有奇异值都很小, 说明函数非常接近 0
    alpha = [];
    c = [];
    return;
end

% 4. 构造 U^(+) 和 U^(-)，并解特征值问题
% --------------------------------------
% 令 U^(+)=U_r(1:N,:), U^(-)=U_r(2:N+1,:)
Uplus = U_r(1:N, :);      % 大小: N x r
Uminus = U_r(2:N+1, :);   % 大小: N x r

% 矩阵 pencil 方法: Z = pinv(Uplus) * Uminus
Z = pinv(Uplus) * Uminus; % 尺寸: r x r

% 解特征值问题
[EigVec, EigValMat] = eig(Z);  % Z*w = lambda*w
lambda = diag(EigValMat);     

dt = T / (2*N);
alpha = (1/dt) * log(lambda);
Vmat = exp(tSamples_pos*(alpha.'));

% 求解 c, 使得 Vmat * c ~ ySamples
c = pinv(Vmat)*ySamples_pos; 
alpha = double(alpha);
c = double(c);

%%
d = 600;
s = -alpha;
w = c;
s = mp(s,d);
A = diag(-s);
B = sqrt(mp(w,d));
C = B;
C = C.';
n = length(s);
T = 10;
P = B*B';
Q = C'*C;
AA = mp(s + s',d);
AAA = AA.^2;
EA = mp(exp(-T*s),d);
%  w(r) = 1
Ptl = P.*(1- EA*EA')./(AA);
Qtl = conj(Ptl);

% w(r) = 1./sqrt(r + K)
K = 1e-2;
EWA = exp(K*AA);
ggt = expint(AA*(T+K));
gg0 = expint(AA*(K));
P = P.*EWA.*(gg0-ggt);
Q = conj(P);

% weight = self
% EP = mp(zeros(length(s),length(s)));
% for i = 1:length(s)
%     AAE = mp(s + s'+s(i),d);
%     EAE = mp(exp(-T*(s+s(i)/2)),d);
%     EP = EP + w(i)*(1- EAE*EAE')./(AAE);
% end
% P = P.*EP;
% Q = conj(P);


%%
% for i=1:n-1
%     P(i,i)=sqrt(P(i,i)-sum(P(i,1:i-1).^2));
%     Q(i,i)=sqrt(Q(i,i)-sum(Q(i,1:i-1).^2));
%     PP=sum(P(i+1:n-1,1:i-1).*P(i,1:i-1),2);
%     QQ=sum(Q(i+1:n-1,1:i-1).*Q(i,1:i-1),2);
%     P(i+1:n-1,i)=(P(i+1:n-1,i)-PP(1:n-1-i))/P(i,i);
%     Q(i+1:n-1,i)=(Q(i+1:n-1,i)-QQ(1:n-1-i))/Q(i,i);
% end
% P=tril(ones(n-1)).*P;Q=tril(ones(n-1)).*Q;
% Lc=P;
% LL=P.'*Q;
RPtl = chol(Ptl,'lower');
RQtl = chol(Qtl,'lower');

RP = chol(P,'lower');
RQ = chol(Q,'lower');




Stl = RPtl;
Ltl = RQtl;
S = RP;
L = RQ;
%% 
LLtl = Stl'*Ltl;
LL = S'*L;
[Utl,Sigmatl,~] = svd(LLtl);
[U,Sigma,~] = svd(LL);

%% 
LLL=mp(diag(mp(Sigmatl,d)),d);
LLL=diag(mp(LLL.^(mp(-1/2,d)),d));
Ttl = mp(Stl*Utl*LLL,d);
invTtl = inv(Ttl);

LLL=mp(diag(mp(Sigma,d)),d);
LLL=diag(mp(LLL.^(mp(-1/2,d)),d));
T = mp(S*U*LLL,d);
invT = inv(T);
%% 
Atl = mp(invTtl,d)*mp(A,d)*mp(Ttl,d);
Btl = mp(invTtl,d)*mp(B,d);
Ctl = mp(C,d)*mp(Ttl,d);

At = mp(invT,d)*mp(A,d)*mp(T,d);
Bt = mp(invT,d)*mp(B,d);
Ct = mp(C,d)*mp(T,d);

%% 
T = 10;
tGrid =linspace(0, T,1000000);
fExact = fTest(tGrid); 
maxerror_list_tlbt = n1:r;
maxerror_list_wbt = n1:r;

%% 
p=27;
Adtl = Atl(1:p, 1:p);
Bdtl = Btl(1:p);
Cdtl = Ctl(1:p);


Ad = At(1:p, 1:p);
Bd = Bt(1:p);
Cd = Ct(1:p);

[Vtl,Dtl,W] = eig(Adtl);
[V,D,W] = eig(Ad);
s_tl = - diag(Dtl);
s_nr = - diag(D);
s_tl = s_tl.';
s_nr = s_nr.'; % new s

B_tl = inv(Vtl)*Bdtl;
C_tl = Cdtl*Vtl;
w_tl = B_tl.*(C_tl.');

B_mr = inv(V)*Bd;
C_mr=Cd*V;
w_nr = B_mr.*(C_mr.');

s_nr = double(s_nr);
w_nr = double(w_nr);

s_tl = double(s_tl);
w_tl = double(w_tl);
% save(filename, "s_d","w_d");  % save as double
fApprox_tlbt = zeros(size(tGrid));
for k = 1:p
    fApprox_tlbt = fApprox_tlbt + w_tl(k) * exp(-s_tl(k) * tGrid);
end

fApprox_wbt = zeros(size(tGrid));
for k = 1:p
    fApprox_wbt = fApprox_wbt + w_nr(k) * exp(-s_nr(k) * tGrid);
end
absError_tlbt = abs(fExact - real(fApprox_tlbt));
absError_wbt = abs(fExact - real(fApprox_wbt));

plot(tGrid,log10(absError_wbt));
hold on;
plot(tGrid,log10(absError));
plot(tGrid,log10(absError_tlbt));
xlim([0,10])
xlabel('p')
ylabel('log10(Maximum AbsError)')
legend('WBT', 'HSVD', 'TLBT', 'Location', 'Best');  % 'Location' 参数可以调整图例在图中的位置
hold off;
hold off
%% 

for p=n1:r
    Adtl = Atl(1:p, 1:p);
    Bdtl = Btl(1:p);
    Cdtl = Ctl(1:p);


    Ad = At(1:p, 1:p);
    Bd = Bt(1:p);
    Cd = Ct(1:p);
    
    [Vtl,Dtl,W] = eig(Adtl);
    [V,D,W] = eig(Ad);
    s_tl = - diag(Dtl);
    s_nr = - diag(D);
    s_tl = s_tl.';
    s_nr = s_nr.'; % new s

    B_tl = inv(Vtl)*Bdtl;
    C_tl = Cdtl*Vtl;
    w_tl = B_tl.*(C_tl.');

    B_mr = inv(V)*Bd;
    C_mr=Cd*V;
    w_nr = B_mr.*(C_mr.');
    
    s_nr = double(s_nr);
    w_nr = double(w_nr);

    s_tl = double(s_tl);
    w_tl = double(w_tl);
    % save(filename, "s_d","w_d");  % save as double
    fApprox_tlbt = zeros(size(tGrid));
    for k = 1:p
        fApprox_tlbt = fApprox_tlbt + w_tl(k) * exp(-s_tl(k) * tGrid);
    end

    fApprox_wbt = zeros(size(tGrid));
    for k = 1:p
        fApprox_wbt = fApprox_wbt + w_nr(k) * exp(-s_nr(k) * tGrid);
    end
    
    absError_tlbt = abs(fExact - real(fApprox_tlbt));
    absError_wbt = abs(fExact - real(fApprox_wbt));
    maxerror_list_tlbt(p-n1+1) = max(absError_tlbt);
    maxerror_list_wbt(p-n1+1) = max(absError_wbt);
    % plot(tGrid,log10(absError))
    % disp(max(absError))
end
%% 

plot([(12):2:r,r],log10(maxerror_list([(12):2:r,r])), 'r-', 'LineWidth', 2);
hold on;
plot([(12):2:r,r],log10(maxerror_list_tlbt([(12):2:r,r])), 'b-', 'LineWidth', 2);
plot([(12):2:r,r],log10(maxerror_list_wbt([(12):2:r,r])), 'g-', 'LineWidth', 2);
% title(sprintf('HSVD vs HSVD-WBT'))
xlim([12,r])
xlabel('p')
ylabel('log10(Maximum AbsError)')
legend('HSVD', 'HSVD+TLBT', 'HSVD+WBT', 'Location', 'Best');  % 'Location' 参数可以调整图例在图中的位置
hold off;



function [alpha, c, r] = HSVD(fHandle,shift, N, T, r)
%  使用 Hankel 矩阵 + SVD + 小规模特征值问题，得到对函数 f(t)
%  的有限个指数和近似表示：
%      f(t) ~ sum_{k=1}^r c_k * exp(alpha_k * t).
%
% 输入:
%   fHandle : 函数句柄, 例如 @(t) exp(-t) + sin(t)
%   N       : 一半的采样数量 - 1, 整体采样点数 = 2N+1
%   T       : 采样区间 [0, T]
%   r     : 截断项数, 用于确定缩减后指数项数目
%
% 输出:
%   alpha   : 大小为 r 的列向量, 存储每个指数项 e^{alpha_k t} 中的 alpha_k
%   c       : 大小为 r 的列向量, 存储每个指数项的系数 c_k
%   r       : 实际保留的指数项数目(截断后秩)
%
% 例: [alpha, c, r] = computeExponentialRepresentation(@(t) exp(-2*t)+0.1*sin(5*t), 50, 10, 1e-10);

    % 1. 采样
    % --------
    % 在 [0, T] 区间等距采样 2N+1 个点
    tSamples = linspace(-shift, T, 2*N + 1).';
    tSamples_pos = tSamples(tSamples >= 0);
    % ySamples = fHandle(tSamples);
    ySamples_pos = fHandle(tSamples_pos);
    H = tSamples + tSamples.';
    H = H(1:N+1,1:N+1);
    H = fHandle(H);
    
    % 3. 对 Hankel 矩阵做 SVD, 并根据 tol 截断
    % --------------------------------------
    [U, S, V] = svd(H, 'econ');      % H = U*S*V',  'econ'可节省空间
    sVals = diag(S);                 % 奇异值向量
    sMax = sVals(1);                 % 最大奇异值
    % idx = find(sVals > tol * sMax);  % 保留的奇异值下标
    % r = length(idx);                 % 截断后实际秩
    U_r = U(:, 1:r);
    
    if r == 0
        % 如果所有奇异值都很小, 说明函数非常接近 0
        alpha = [];
        c = [];
        return;
    end
    
    % 4. 构造 U^(+) 和 U^(-)，并解特征值问题
    % --------------------------------------
    % 令 U^(+)=U_r(1:N,:), U^(-)=U_r(2:N+1,:)
    Uplus = U_r(1:N, :);      % 大小: N x r
    Uminus = U_r(2:N+1, :);   % 大小: N x r
    
    % 矩阵 pencil 方法: Z = pinv(Uplus) * Uminus
    Z = pinv(Uplus) * Uminus; % 尺寸: r x r
    
    % 解特征值问题
    [EigVec, EigValMat] = eig(Z);  % Z*w = lambda*w
    lambda = diag(EigValMat);     
    
    dt = T / (2*N);
    alpha = (1/dt) * log(lambda);
    Vmat = exp(tSamples_pos*(alpha.'));
    
    % 求解 c, 使得 Vmat * c ~ ySamples
    c = pinv(Vmat)*ySamples_pos; 
end

