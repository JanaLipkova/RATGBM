%========================================
%
%   Correlation plot between proliferation rate from:
%   - UQ inference
%   - histology inference
%
%========================================



%plotting set up
fs=30;
lw = 2;

rhoRange = [2.7e-03, 1.9e-01]; %[1/day]

% Rho from UQ inference Mean /pm std
P01rhoUQ = [ 0.0233    0.0845];
P07rhoUQ = [ 0.0027    0.0693];
P11rhoUQ = [ 0.0027    0.0370];
P22rhoUQ = [ 0.0038    0.0466];

% MAP \pm std
% P01rhoUQ = [ max(0.0288 - 0.0306, rhoRange(1)),   min(0.0288 + 0.0306, rhoRange(2)) ];
% P07rhoUQ = [ max(0.0238 - 0.0336, rhoRange(1)),   min(0.0238 + 0.0336, rhoRange(2)) ];
% P11rhoUQ = [ max(0.0133 - 0.0174, rhoRange(1)),   min(0.0133 + 0.0174, rhoRange(2)) ];
% P22rhoUQ = [ max(0.00846 - 0.02144, rhoRange(1)),   min(0.00846 + 0.02144, rhoRange(2)) ];



% 95 % CI
% P01rhoUQ = [ 0.0173, 0.1388];
% P07rhoUQ = [ 0.0084   0.1396];
% P11rhoUQ = [0.0072,  0.0851];
% P22rhoUQ = [ 0.0059,  0.0879];

%65 % CI
% P01rhoUQ = [0.0280,0.0770];
% P07rhoUQ = [ 0.0119   0.0588];
% P11rhoUQ = [0.0109, 0.0229];
% P22rhoUQ = [ 0.0101,  0.03719];


% Rho from Histology (assume cell cycle of 4-8 days) 
P01rhoH = [ 0.0228     0.0456];
P07rhoH = [ 0.0037     0.0074];
P11rhoH = [ 0.0248     0.0558];
P22rhoH = [ 0.0248     0.0558];

% Points of the squres
P01x = [ P01rhoUQ(1), P01rhoUQ(2), P01rhoUQ(2), P01rhoUQ(1), P01rhoUQ(1)];
P01y = [ P01rhoH(1),  P01rhoH(1),  P01rhoH(2),  P01rhoH(2),  P01rhoH(1) ];

P07x = [ P07rhoUQ(1), P07rhoUQ(2), P07rhoUQ(2), P07rhoUQ(1), P07rhoUQ(1)];
P07y = [ P07rhoH(1),  P07rhoH(1),  P07rhoH(2),  P07rhoH(2),  P07rhoH(1) ];

P11x = [ P11rhoUQ(1), P11rhoUQ(2), P11rhoUQ(2), P11rhoUQ(1), P11rhoUQ(1)];
P11y = [ P11rhoH(1),  P11rhoH(1),  P11rhoH(2),  P11rhoH(2),  P11rhoH(1) ];

P22x = [ P22rhoUQ(1), P22rhoUQ(2), P22rhoUQ(2), P22rhoUQ(1), P22rhoUQ(1)];
P22y = [ P22rhoH(1),  P22rhoH(1),  P22rhoH(2),  P22rhoH(2),  P22rhoH(1) ];

% Compute mid points for correlations
P01midUQ = P01rhoUQ(1) + 0.5 * (P01rhoUQ(2) - P01rhoUQ(1));
P07midUQ = P07rhoUQ(1) + 0.5 * (P07rhoUQ(2) - P07rhoUQ(1));
P11midUQ = P11rhoUQ(1) + 0.5 * (P11rhoUQ(2) - P11rhoUQ(1));
P22midUQ = P22rhoUQ(1) + 0.5 * (P22rhoUQ(2) - P22rhoUQ(1));

P01midH = P01rhoH(1) + 0.5 * (P01rhoH(2) - P01rhoH(1));
P07midH = P07rhoH(1) + 0.5 * (P07rhoH(2) - P07rhoH(1));
P11midH = P11rhoH(1) + 0.5 * (P11rhoH(2) - P11rhoH(1));
P22midH = P22rhoH(1) + 0.5 * (P22rhoH(2) - P22rhoH(1));




figure;
set(gca,'Fontsize',fs);
hold on

rhoRange = rhoRange./2.2;
 axis([rhoRange(1),rhoRange(2),rhoRange(1), rhoRange(2)])
% set(gca,'XScale','log');
% set(gca,'YScale','log');
set(gca,'linewidth',3) 
axis square;

plot(P01x, P01y, 'b--', 'LineWidth', 3);
plot(P07x, P07y, 'g--', 'LineWidth', 3);
plot(P11x, P11y, 'k--', 'LineWidth', 3);
plot(P22x, P22y, 'y--', 'LineWidth', 3);


% legend('P01','P07','P11','P22')
legend('P01','P02','P03','P04')
plot( rhoRange,rhoRange,'r-','LineWidth', 3)

title('Mean \pm std ')

xlabel('\rho^{UQ} [1/day]');
ylabel('\rho^{H} [1/day]');

box on;  grid on

set(gcf, 'Position', [100, 100, 1049, 895]);
set(gcf,'papersize',[400,400]);
set(gcf,'PaperPositionMode','auto')
print('Histo_MeanSTD_Linear','-djpeg')


