%%% plot_lab_results

krytar = [3.45 4.48 5.52 6.52 7.67 8.67];
correction = 8.05 - 8.70;
ber = [1.9e-3 7.0e-4 1.7e-4 3.5e-5 6e-6 1.3e-6];

% x = [-65.42 -66.25 -67.46 -68.06];
% y = [5.5e-5 4.37e-5 7.5e-3 4.13e-2];

x = krytar + correction - 80;
y = ber;

figure(99);
semilogy(x,y,'.-'); grid on;
xlabel('input signal level (dBm)');
ylabel('BER');

%%% Matlab Version < 2015a

% homer = get(gcf,'Units');
% set(gcf,'Units','inches');
% bart = get(gcf,'Position');
% set(gcf,'Position',[bart(1:2) 5 4]);
% set(gcf,'PaperPositionMode','auto');
% set(gcf,'Units',homer);

%%% Matlab Version >= 2015a

f = gcf;
homer = f.Units;
f.Units = 'inches';
bart = f.Position;
f.Position = [bart(1:2) 5 4];
f.PaperPositionMode = 'auto';
f.Units = homer;

print -depsc lab-test.eps