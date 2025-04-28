clc
clear
format long
% d =1000;
% mp(1000);

eta = [2.1073876854180e-12, 1.8365780986634e-11, 4.7777245228151e-11, 8.5624630300630e-11, 1.3289239111902e-10, 2.0054640049463e-10, 3.0217586807074e-10, 4.5529860118663e-10];
weta = [3.2630674210379e-6, 3.1058837221013e-6, 2.8014247111005e-6, 2.5227064974618e-6, 2.7039982943831e-6, 3.2761422288967e-6, 4.0205002817225e-6, 4.9351231646262e-6];
% eta = mp(eta,d);
% w = mp(w,d);
s_wbt = [4.551547331769476e-10,2.967225833661697e-10,1.69007319959171e-10,7.549432548035814e-12,6.564737578696118e-11];
w_wbt = [4.955235574250308e-6,4.503807233967136e-6, 5.145266724826999e-6, 6.175199783823309e-6, 5.849337004446485e-6];


h = 0.40994422603935795;
tau = 0.192967891816239;
% h = mp(h,d);
% tau = mp(tau,d);

num = -203:86;
% num = mp(num,d);
% w_rest = h./gamma(1/2).*mp(exp((h.*num-tau)/2),d);
% s_rest = mp(exp(h.*num - tau),d);
w_rest = h./gamma(1/2).*exp((h.*num-tau)/2);
s_rest = exp(h.*num - tau);
% s = [eta,s_rest];
% w = [weta,w_rest];
% s_wbt = [s_wbt,s_rest];
% w_wbt = [w_wbt, w_rest];
s = s_rest;
w = w_rest;
%% 
x = [];
for i=-7:3
    x = [x,linspace(10^(i),10^(i+1),1000)];
end
x = [x,linspace(10^(4),10^(5),10000)];
%% 
x = x.';
xx = x.^2;
%% 

sog = exp(-xx*s)*w.';
% sog_wbt = exp(-xx*s_wbt)*w_wbt.';
error = abs(sog-1./x);
% error_wbt = abs(sog_wbt - 1./x);

%% 

rel_beylkin_error = error.*x;
loglog(x,rel_beylkin_error,'r-', 'LineWidth', 2);

%% MR from Beylkin
d = 500;
s = s_rest(1:47);
w = w_rest(1:47);
sog0 = exp(-xx*s)*w.';
loglog(x,sog0);
%% 
s = s;
ww = w;
A = diag(-s);
B = sqrt(mp(abs(ww),d)).';
n = length(s);
% T = (1e10)/2;
P = B*B.';
AA = mp(s.' + s,d);
AAA = AA.^2;
% EA = mp(exp(-T*s),d);
% BT
P = P.*(1)./AA;


% %  w(r) = 1
% P = P.*(1 - EA.'*EA)./AA;

% % w(r) = sqrt(r)
% P = P.*(1 - (EA.'*EA).*(AA*T + 1))./AAA;

% w(r) = r
% AAAA = AA.^3;
% P = P.*(2 - (EA.'*EA).*(AA*T.*(AA*T+2) + 2))./AAAA;

% % w(r) = sqrt(r + 1)
% P = P.*(1 + AA - (EA.'*EA).*(AA*T + AA + 1))./AAA;

% w(r) = 1./sqrt(r + K)
% K = 1e-4;
% EWA = exp(K*AA);
% ggt = expint(AA*(T+K));
% gg0 = expint(AA*(K));
% P = P.*EWA.*(gg0-ggt);

% 分段
% K = 1e+8;
% EKA = mp(exp(-K*s),d);
% P1 = 1e-8*P.*(1 - (EKA.'*EKA).*(AA*K + 1))./AAA;
% EWA = exp(K*AA);
% ETKA = mp(exp(-(T-K)*s),d);
% ggt = expint(AA*(T));
% gg0 = expint(AA*(K));
% P2 = K*P.*(ETKA.'*ETKA) .*(EWA.*(gg0-ggt));
% P = (P1+P2);

%% 
[U,S,V] = svd(P);

At = U.'*A*U;
Bt = U'* B;

%% Rebuild
p = 46;
Ad = At(1:p, 1:p);
Bd = Bt(1:p);

[V,D,W] = eig(Ad);

s_nr = - diag(D);
s_nr = s_nr.'; % new s
s_nr = s_nr;
B_mr = Bd'*V;

w_nr = B_mr.^2;
w_nr = w_nr.'; % new w
% w_nr = w_nr.*mp(exp(mp(1e+1*s_nr.',d)),d);
% filename = sprintf('%s/p_%d_mp.mat', folder, p);
% save(filename, 's_nr', 'w_nr'); % save as high precision

% filename = sprintf('%s/p_%d_double.mat', folder, p);
% s_d = double(s_nr);
% w_d = double(w_nr);
% save(filename, "s_d","w_d");  % save as double

sog1 = exp(mp(-xx*s_nr,d))*w_nr;
% sog0 = exp(mp(-xx*s,d))*w.';
error = abs(sog1- sog0);% alpha = 1
rel_error = error.*x;

loglog(x,rel_error)
xlabel("r")
ylabel("Relative Error")
%% 

loglog(xx,sog0)

%% 
s_new = [s_rest(153:end),s_nr];
w_new = [w_rest(153:end).';w_nr];

s_new = double(s_new);
w_new = double(w_new);

sog_new = exp(-xx*s_new)*w_new;
error = abs(sog_new - 1./x);
rel_error_total = error.*x;

loglog(x,rel_error_total)

% n = 20;
% r_error = zeros(1,20);
% folder = 'nearly_singular_WBT';
% if ~exist(folder, 'dir')
%     mkdir(folder);
% end
% for p = 81:100
%     Ad = At(1:p, 1:p);
%     Bd = Bt(1:p);
% 
%     [V,D,W] = eig(Ad);
% 
%     s_nr = - diag(D);
%     s_nr = s_nr.'; % new s
%     B_mr = Bd'*V;
% 
%     w_nr = B_mr.^2;
%     w_nr = w_nr.'; % new w
%     filename = sprintf('%s/p_%d_mp.mat', folder, p);
%     save(filename, 's_nr', 'w_nr'); % save as high precision
% 
%     filename = sprintf('%s/p_%d_double.mat', folder, p);
%     s_d = double(s_nr);
%     w_d = double(w_nr);
%     save(filename, "s_d","w_d");  % save as double
% 
%     sog1 = mp(exp(mp(-xx*s_nr,d)),d)*w_nr;
%     merror = abs(sog1- 1 ./(x));% alpha = 1
%     rel_error = merror.*x;
%     r_error(p-80) = max(rel_error);% save maximum error
% end

%% 
figure;

% 绘制第一条 loglog 曲线，设置线型和线宽
loglog(x, rel_error_1, 'r-', 'LineWidth', 2);  % 红色实线

% 使用 hold on 保持当前图形，方便继续在同一图内绘图
hold on;

% 绘制第二条曲线
% loglog(x, rel_error_sqrtr, 'b--', 'LineWidth', 2); % 蓝色虚线

% 绘制第三条曲线
loglog(x, rel_error_1sqrtre7, 'g-.', 'LineWidth', 2); % 绿色点划线

% 绘制第三条曲线
loglog(x, rel_error_1sqrtre10, 'p-.', 'LineWidth', 2); % 绿色点划线

% 绘制第三条曲线
loglog(x, rel_error_exp10, 'k--', 'LineWidth', 2); % 绿色点划线

% 添加轴标签、标题以及图例
xlabel('X轴');
ylabel('Relative Error');
legend('w(r)=1', 'w(r)=1/sqrt(r+1e-7)', 'w(r)=1/sqrt(r+1e-10)','w(r)=exp(-5r)', 'Location', 'Best');  % 'Location' 参数可以调整图例在图中的位置

% 如果不再绘制其他内容，记得关闭 hold 状态
hold off;

