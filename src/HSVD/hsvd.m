format long
global d
T = 100;         % 计算区间 [0, T]
N = 8000;         % 采样参数 (总采样点为 2N+1)
% tol = 1e-10;    % 奇异值截断阈值
r = 40;
L = 100;
shift = 1;
d = 200;
mp(d);

fHandle = @(t) erf(L*(t+1e-16))./(t+1e-16);  
fTest = @(t) erf(L*(t+1e-16))./(t+1e-16);
% fHandle = @(t) exp(-t.^2);  
% fTest = @(t) exp(-t.^2);
%% Terms double

n1 = 10; n2 = 25;
maxerror_list = n1:n2;
for r= n1:n2
    [alpha, c, r] = HSVD(fTest,shift, N, T, r);
    tGrid =linspace(0, T,10000);    % 在 [0,T] 上取 1000 个点进行评估
    fExact = fTest(tGrid);           % 真值
    alpha = double(alpha);
    c = double(c);
    % 计算近似值: fApprox(t) = sum_{k=1}^{r} c_k * exp(alpha_k*t)
    fApprox = mp(zeros(size(tGrid)),d);
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
plot(n1:n2,log10(maxerror_list));
title(sprintf('Lambda = %.3f, shift = %s, Sample num = %.3f',L , string(shift), 2*N+1))
xlabel('p')
ylabel('log10(Maximum AbsError)')
%% 
plot(n1:n2,log10(maxerror_list05), 'r-', 'LineWidth', 2);
hold on;
plot(n1:n2,log10(maxerror_list1), 'g-', 'LineWidth', 2);
plot(n1:n2,log10(maxerror_list2), 'b-', 'LineWidth', 2);
plot(n1:n2,log10(maxerror_list4), 'k-', 'LineWidth', 2);
title(sprintf('HSVD on Ewald kernel with %.f sample points', 2*N+1))
xlabel('p')
ylabel('log10(Maximum AbsError)')
legend('\Lambda = 0.5', '\Lambda = 1', '\Lambda = 2','\Lambda = 4', 'Location', 'Best');  % 'Location' 参数可以调整图例在图中的位置
hold off;

%% 
r=28;
[alpha, c, r] = HSVD(fTest,shift, N, T, r); 

tGrid =linspace(0, T,100000);    % 在 [0,T] 上取 1000 个点进行评估
fExact = fTest(tGrid);           % 真值
alpha = double(alpha);
c = double(c);
% 计算近似值: fApprox(t) = sum_{k=1}^{r} c_k * exp(alpha_k*t)
fApprox = mp(zeros(size(tGrid)),d);
for k = 1:r
    fApprox = fApprox + c(k) * exp(alpha(k) * tGrid);
end

absError = abs(fExact - real(fApprox));
relError = abs(fExact - real(fApprox)) ./ abs(fExact);
maxRelError = max(relError);
maxAbsError = max(absError);
plot(tGrid,absError);
title(sprintf('Lambda = %.1f, p = %.1f, shift = %s',L ,r, string(shift)))
xlabel('r')
ylabel('absError')

%% 
s = -alpha;
w = c;
A = diag(-s);
B = sqrt(mp(abs(w),d));
C = mp(sign(w),d).*B;
C = C.';
n = length(s);
T = 10;
P = B*B';
Q = C'*C;
AA = mp(s' + s,d);
AAA = AA.^2;
EA = mp(exp(-T*s),d);
%  w(r) = 1
P = P.*(1 - (EA'*EA).')./(AA.');
Q = Q.*(1 - EA'*EA)./AA;
%% 
n=17;
for i=1:n-1
    P(i,i)=sqrt(P(i,i)-sum(P(i,1:i-1).^2));
    Q(i,i)=sqrt(Q(i,i)-sum(Q(i,1:i-1).^2));
    PP=sum(P(i+1:n-1,1:i-1).*P(i,1:i-1),2);
    QQ=sum(Q(i+1:n-1,1:i-1).*Q(i,1:i-1),2);
    P(i+1:n-1,i)=(P(i+1:n-1,i)-PP(1:n-1-i))/P(i,i);
    Q(i+1:n-1,i)=(Q(i+1:n-1,i)-QQ(1:n-1-i))/Q(i,i);
end
P=tril(ones(n-1)).*P;Q=tril(ones(n-1)).*Q;
Lc=P;
LL=P.'*Q;
%% 


[U,S,V] = svd(LL);
LLL=mp(diag(mp(S,d)),d);
LLL=mp(LLL.^(mp(-1/2,d)),d);
T=mp(mp(Lc,d)*mp(U,d)*mp(diag(LLL),d),d);
% InvT = mp(diag(LLL),d)*mp(V',d)*mp(Q',d);
InvT = inv(T);
I = T*T.';
At = mp(InvT,d)*mp(A,d)*mp(T,d);
Bt = mp(InvT,d)*mp(B,d);
Ct = mp(C,d)*mp(T,d);
%% 

x = tGrid.';
p = 15;
Ad = At(1:p, 1:p);
Bd = Bt(1:p);
Cd = Ct(1:p);

[V,D,W] = eig(Ad);
%% 

s_nr = - diag(D);
s_nr = s_nr.'; % new s
B_mr = W*Bd;
C_mr=Cd*V;



w_nr = B_mr.*(C_mr.');
% w_nr = w_nr.'; % new w
% filename = sprintf('%s/p_%d_mp.mat', folder, p);
% save(filename, 's_nr', 'w_nr'); % save as high precision
% 
% filename = sprintf('%s/p_%d_double.mat', folder, p);
s_nr = double(s_nr);
w_nr = double(w_nr);
% save(filename, "s_d","w_d");  % save as double
%% 

fApprox_wbt = zeros(size(tGrid));
for k = 1:p
    fApprox_wbt = fApprox_wbt + w_nr(k) * exp(-s_nr(k) * tGrid);
end

absError = abs(fExact - real(fApprox_wbt));
plot(tGrid,log10(absError))

%% Single MP

% [alpha, c, r] = HSVD(fTest,shift, N, T, r);
% tGrid =linspace(0, T,100000);    % 在 [0,T] 上取 1000 个点进行评估
% fExact = fTest(tGrid);           % 真值
% alpha = double(alpha);
% c = double(c);
% % 计算近似值: fApprox(t) = sum_{k=1}^{r} c_k * exp(alpha_k*t)
% fApprox = mp(zeros(size(tGrid)),d);
% for k = 1:r
%     fApprox = fApprox + c(k) * exp(alpha(k) * tGrid);
% end
% 
% absError = abs(fExact - real(fApprox));
% relError = abs(fExact - real(fApprox)) ./ abs(fExact);
% 
% % 计算在整个区间上的最大相对误差
% maxRelError = max(relError);
% maxAbsError = max(absError);
% plot(tGrid,log10(absError));
% title(sprintf('Lambda = %.3f, p = %.1f, shift = %s, Sample num = %.3f',L ,r, string(shift), 2*N+1))
% xlabel('x')
% ylabel('log10(Maximum AbsError)')

% [alpha, c, r] = HSVD_mp(fTest,shift, N, T, r);
% tGrid =linspace(0, T,10000);    % 在 [0,T] 上取 1000 个点进行评估
% fExact = fTest(tGrid);           % 真值
% alpha = double(alpha);
% c = double(c);
% % 计算近似值: fApprox(t) = sum_{k=1}^{r} c_k * exp(alpha_k*t)
% fApprox = mp(zeros(size(tGrid)),d);
% for k = 1:r
%     fApprox = fApprox + c(k) * exp(alpha(k) * tGrid);
% end
% 
% absError = abs(fExact - real(fApprox));
% relError = abs(fExact - real(fApprox)) ./ abs(fExact);
% 
% % 计算在整个区间上的最大相对误差
% maxRelError = max(relError);
% maxAbsError = max(absError);
% plot(tGrid,absError);
% title(sprintf('Lambda = %.1f, p = %.1f, shift = %s',L ,r, string(shift)))
% xlabel('r')
% ylabel('AbsError')
% 
%% 

% 
% tSamples = mp(linspace(0, T, 2*N + 1),d).'; 
% ySamples = fHandle(tSamples);
% H = tSamples + tSamples.';
% H = H(1:N+1,1:N+1);
% H = fHandle(H);
% 
% 
% %% 
% 
% [U, S, V] = svd(H, 'econ');     
% sVals = diag(S);                
% sMax = sVals(1);                 % 最大奇异值
% % idx = find(sVals > tol * sMax);  % 保留的奇异值下标
% % r = length(idx);                 % 截断后实际秩
% U_r = U(:, 1:r);  % (N+1) x r
% 
% Uplus = U_r(1:N, :);     
% Uminus = U_r(2:N+1, :);   
% 
% Z = pinv(Uplus) * Uminus;
% 
% [EigVec, EigValMat] = eig(Z);  % Z*w = lambda*w
% lambda = diag(EigValMat);     
% 
% dt = mp(T / (2*N),d);
% alpha = (1/dt) * mp(log(lambda),d);
% Vmat = mp(exp(tSamples*(alpha.')),d);
% 
% % 求解 c, 使得 Vmat * c ~ ySamples
% c = pinv(Vmat)*ySamples; 
% %% 

% disp('指数项数目 r = ');
% disp(r);
% 
% disp('alpha = ');
% disp(alpha);
% 
% disp('c = ');
% disp(c);


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
    ySamples = fHandle(tSamples);
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
    Vmat = exp(tSamples*(alpha.'));
    
    % 求解 c, 使得 Vmat * c ~ ySamples
    c = pinv(Vmat)*ySamples; 
end


function [alpha, c, r] = HSVD_mp(fHandle,shift, N, T, r)
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
    global d
    tSamples = mp(linspace(-shift, T, 2*N + 1),d).'; 
    ySamples = fHandle(tSamples);
    H = tSamples + tSamples.';
    H = H(1:N+1,1:N+1);
    H = fHandle(H);
    
    [U, S, V] = svd(H, 'econ');     
    sVals = diag(S);                
    sMax = sVals(1);                 % 最大奇异值
    % idx = find(sVals > tol * sMax);  % 保留的奇异值下标
    % r = length(idx);                 % 截断后实际秩
    U_r = U(:, 1:r);  % (N+1) x r
    
    Uplus = U_r(1:N, :);     
    Uminus = U_r(2:N+1, :);   
    
    Z = pinv(Uplus) * Uminus;
    
    [EigVec, EigValMat] = eig(Z);  % Z*w = lambda*w
    lambda = diag(EigValMat);     
    
    dt = mp(T / (2*N),d);
    alpha = (1/dt) * mp(log(lambda),d);
    Vmat = mp(exp(tSamples*(alpha.')),d);
    
    % 求解 c, 使得 Vmat * c ~ ySamples
    c = pinv(Vmat)*ySamples; 
    
end