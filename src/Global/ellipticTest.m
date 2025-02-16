clear all;
x = -0.5:0.05:1.5;
m = sqrt(1-x);
n = length(x);
k = zeros(n,1);
e = zeros(n,1);
p = zeros(n,1);
km = zeros(n,1);
em = zeros(n,1);
pm = zeros(n,1);

for i=1:n
    [k(i),e(i),p(i)] = elliptic123(x(i),0.2);
    km(i)= ellipticK(x(i));
    em(i) = ellipticE(x(i));
    pm(i) = ellipticPi(0.2, x(i));
    
end

%%
figure;
title('Elliptic Integrals Evaluation');
plot(x,k,'k-','linewidth',1); hold on;
plot(x,e,'b-','linewidth',1); hold on;
plot(x,p,'m-','linewidth',1); hold on;
plot(x,km,'ks','linewidth',1);
plot(x,em,'bs','linewidth',1); hold on;
plot(x,pm,'ms','linewidth',1);
legend('coded K','coded E','coded Pi','Matlab K','Matlab E','Matlab Pi');
set(gca, 'FontName','Times')
set(gca, 'FontSize',23);


