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

num = -51:86;
% num = mp(num,d);
% w_rest = h./gamma(1/2).*mp(exp((h.*num-tau)/2),d);
% s_rest = mp(exp(h.*num - tau),d);
w_rest = h./gamma(1/2).*exp((h.*num-tau)/2);
s_rest = exp(h.*num - tau);
s = [eta,s_rest];
w = [weta,w_rest];
s_wbt = [s_wbt,s_rest];
w_wbt = [w_wbt, w_rest];
% s = s_rest;
% w = w_rest;
%% 
x = [];
for i=-8:3
    x = [x,linspace(10^(i),10^(i+1),1000)];
end
x = [x,linspace(10^(4),10^(5),1000)];
x = [x,linspace(10^(5),10^(6),1000)];
%% 
x = x.';
xx = x.^2;
%% 

sog = exp(-xx*s)*w.';
sog_wbt = exp(-xx*s_wbt)*w_wbt.';
error = abs(sog-1./x);
error_wbt = abs(sog_wbt - 1./x);

%% 

rel_beylkin_error = error.*x;
rel_wbt_error = error_wbt.*x;
% loglog(x,rel_beylkin_error,'b-', 'LineWidth', 2);
% hold on;
loglog(x, rel_wbt_error,'r-', 'LineWidth', 2);
xlabel('r');
ylabel('Relative Error');
xlim([10^(-8),10^(6)]);
% legend('Benchmark', 'WBT', 'Location', 'Best');  % 'Location' 参数可以调整图例在图中的位置
hold off;
