function [s_wbt, w_wbt, error] = WBT(s, w, x, p,T, method, opt)
% WBT  Weighted/Time-limited/classical balanced truncation with general α-weight
%   [s_wbt, w_wbt, error] = WBT(s,w,p,T,x,method,opt)
%   s,w         original poles and weights for SOE: sum(w.*exp(-s*x))
%   p           truncation order (integer)
%   T           integration upper limit
%   x           vector of test points for error computation
%   method      "classical" | "TLBT" | "WBT"
%   opt.alpha   exponent α for weight (r+K)^(-α) (default 0.5)
%   opt.K       shift K for weight (default 10)
    %% 输入声明 & 默认值
    arguments
        s       (:,1) double                   % 原始 s 向量, 其中原始SOE为 Σw*exp(-s*x) 
        w       (:,1) double                   % 原始 w 向量, 其中原始SOE为 Σw*exp(-s*x)
        x       (1,:) double        % 误差测试点
        p       (1,1) double {mustBeInteger}    % 截断阶数，整数
        T       (1,1) double {mustBePositive}   % 积分上限
        method  {mustBeMember(method,["classical","TLBT","WBT"])} = "TLBT"
        opt.alpha   (1,1) double = 0.5              % 使用w(r)=(r+K)^(-alpha) 时的alpha                                         
        opt.K   (1,1) double = 10                   % 使用w(r)=(r+K)^(-alpha) 时的K
    end

    %% 准备矩阵
    A   = diag(-s);
    B   = sqrt(w);
    C   = B.';
    AA  = s + s';
    EA  = exp(-T*s);
    
    % P 初始构造
    P = B*B';
    Q = C'*C;

    switch method
        case "classical"
            % Classical MR:
            P = P./AA;
            Q = Q./(AA.');

        case "TLBT"
            % Time-limited BT: w(t)=1
            P = P.*(1- EA*EA')./(AA);
            Q = Q.*(1 - (EA*EA').')./(AA.');

        case "WBT"
            % Weighted BT: w(r)=(r+K)^(-alpha)
            K = opt.K;
            alpha = opt.alpha;
            if alpha == 0.5
                EWA = exp(K*AA);
                ggt = expint(AA*(T+K));
                gg0 = expint(AA*(K));
                P = P.*EWA.*(gg0-ggt);
                Q = conj(P);
            else
                a    = 1 - 2*alpha;
                u1   = AA * K;
                u2   = AA * (T + K);
                
                Q1   = gammainc(u1, a, 'upper');
                Q2   = gammainc(u2, a, 'upper');
                
                G1   = Q1 .* gamma(a);
                G2   = Q2 .* gamma(a);
                
                % 注意 sK 处的 exp(sK)  => expm(AA*K) 如果 AA 是矩阵，可用 elementwise
                EWA  = exp(AA * K);
                P    = P .* EWA .* (G1 - G2) .* (AA.^(2*alpha-1));
                Q    = conj(P);
            end

        otherwise
            error("Unknown method “%s”", method)
    end

     if any(isinf(P(:))) && any(isinf(Q(:)))
        error('Detected Inf in P or Q; exiting.');
    end
    %% 
    RP  = chol(P + 1e-15*eye(size(P)),'lower');
    RQ  = chol(Q + 1e-15*eye(size(Q)),'lower');
    S   = RP;
    L   = RQ;
    LL  = S'*L;
    
    
    %% SVD 及降阶
    [U,Sigma,~] = svd(LL);
    Sigma    =diag(Sigma);
    fprintf('截断处相对奇异值 %f\n', Sigma(p)/Sigma(1));
    Sigma_12 =diag(Sigma.^(-1/2));
    Trans    = S*U*Sigma_12;
    invT     = inv(Trans);
    At       = invT*A*Trans;
    Bt       = invT*B;
    Ct       = C*Trans;
    
    Ad   = At(1:p,1:p);
    Bd   = Bt(1:p);
    Cd   = Ct(1:p);

    [V,D] = eig(Ad);
    s_wbt = -diag(D);
    B_mr  = inv(V)*Bd;
    C_mr  = Cd*V;
    w_wbt = B_mr.*(C_mr.');
    x = x.';
    y_original = exp(-x*s.')*w;
    y_wbt     = exp(-x*s_wbt.')*w_wbt;
    error = abs(y_wbt - y_original);
    rerror = error ./(abs(y_original));
    merror = max(error);
    mrerror = max(rerror);
    fprintf('最大误差 %f\n', merror);
    fprintf('最大相对误差 %f\n', merror);
end
