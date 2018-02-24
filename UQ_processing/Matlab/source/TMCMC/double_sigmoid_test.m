%% double logistic sigmoid
N=101;
s2=0.05;

alpha = zeros(1,N);
dif = linspace(-1,1,N);

for i =1:N
      tmp = dif(i)*dif(i)/s2;
%      tmp = 0.8/s2;
    alpha(i) =  0.5 + 0.5 * sign(dif(i)) * (1 - exp(-tmp));
end;

figure(2)
set(gca,'Fontsize',15);

hold on
g=plot(dif,alpha,'-*','Linewidth',4);
xlabel('ui-uc')
ylabel('alpha')
box on; grid on
