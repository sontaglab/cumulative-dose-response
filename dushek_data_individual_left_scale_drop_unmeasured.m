% now drop the one for 0.976562499 which really had no data, so should not be "zero"!

% data from Trendel ... Dushek, Perfect adaptation of CD8+ T cell responses to constant antigen input over a wide range of affinities is overcome by costimulation https://doi.org/10.1126/scisignal.aay9363
% reading file and plotting means (fig 1) as well as individual runs
% (fig 2)
%
% see also dushek_data from some commented-out commands
%
clear

myfontsize = 36;
mylinewidth = 4;

a=csvread('dushek_high.csv');
% flip the rows so that the last row (values at 1 concentration) is now first and the last is values at 2000 concentration:
b = flipud(a);

% adding these lines:
b = b(2:end,:);
noconcentrations=11;

concentrations = b(:,1);
ave1 = mean(b(:,2:4)');
ave2 = mean(b(:,5:7)');
ave3 = mean(b(:,8:10)');
ave4 = mean(b(:,11:13)');
ave5 = mean(b(:,14:16)');
ave6 = mean(b(:,17:19)');
ave7 = mean(b(:,20:22)');
ave8 = mean(b(:,23:25)');
% plot the average of the three "experiments", meaning measuring total cummuative output as a function of time; e.g. ave1 is the average of the 3 measurements at time 1, etc
figure(1)
plot(concentrations,[ave1' ave2' ave3' ave4' ave5' ave6' ave7' ave8'],'linewidth',mylinewidth)
set(gca, 'XScale', 'log')
legend('1 hour','2 hours','3 hours','4 hours','5 hours','6 hours','7 hours','8 hours','FontSize', myfontsize,'location','northwest','Fontweight','bold')
xlim([0 2000])
ylim([0 1]) % there seems to be a small error and one value is 1.0016,
            % which messes up scales
xlabel('input','Fontweight','bold')
ylabel('cDR','Fontweight','bold');
title('cDR, average of 3 experiments','FontSize', myfontsize,'Fontweight','bold');
ax = gca;
ax.XAxis.FontSize = myfontsize;
ax.YAxis.FontSize = myfontsize;
axis square
%
% now do individual:
%
% experiment 0317:
exp11 = b(:,02)';
exp12 = b(:,05)';
exp13 = b(:,08)';
exp14 = b(:,11)';
exp15 = b(:,14)';
exp16 = b(:,17)';
exp17 = b(:,20)';
exp18 = b(:,23)';
figure(2)
%
t = tiledlayout(1, 3, "TileSpacing", "tight");
nexttile
%subplot(1,3,1)
plot(concentrations,[exp11' exp12' exp13' exp14' exp15' exp16' exp17' exp18'],'linewidth',mylinewidth)
set(gca, 'XScale', 'log')
%legend('1 hour','2 hours','3 hours','4 hours','5 hours','6 hours','7 hours','8 hours')
axis square
xlim([0 2000])
ylim([0 1]) % there seems to be a small error and one value is 1.0016,
            % which messes up scales
xlabel('input','FontSize', myfontsize,'Fontweight','bold')
ylabel('cDR','FontSize', myfontsize,'Fontweight','bold');
ax = gca;
ax.XAxis.FontSize = myfontsize;
ax.YAxis.FontSize = myfontsize;
% experiment 0406:
exp21 = b(:,03)';
exp22 = b(:,06)';
exp23 = b(:,09)';
exp24 = b(:,12)';
exp25 = b(:,15)';
exp26 = b(:,18)';
exp27 = b(:,21)';
exp28 = b(:,24)';
%
%figure(3)
nexttile
%subplot(1,3,2)
plot(concentrations,[exp21' exp22' exp23' exp24' exp25' exp26' exp27' exp28'],'linewidth',mylinewidth)
set(gca, 'XScale', 'log')
legend('1 hour','2 hours','3 hours','4 hours','5 hours','6 hours','7 hours','8 hours','FontSize', myfontsize,'location','northwest','Fontweight','bold');
axis square
xlim([0 2000])
ylim([0 1])
xlabel('input','FontSize', myfontsize,'Fontweight','bold')
yticks([])
%ylabel('cDR') % no room if in tight subplot mode
ax = gca;
ax.XAxis.FontSize = myfontsize;
ax.YAxis.FontSize = myfontsize;
% experiment 0325:
exp31 = b(:,04)';
exp32 = b(:,07)';
exp33 = b(:,10)';
exp34 = b(:,13)';
exp35 = b(:,16)';
exp36 = b(:,19)';
exp37 = b(:,22)';
exp38 = b(:,25)';
%
%figure(4)
nexttile
%subplot(1,3,3)
plot(concentrations,[exp31' exp32' exp33' exp34' exp35' exp36' exp37' exp38'],'linewidth',mylinewidth)
set(gca, 'XScale', 'log')
%legend('1 hour','2 hours','3 hours','4 hours','5 hours','6 hours','7 hours','8 hours')
axis square
xlim([0 2000])
ylim([0 1])
xlabel('input','FontSize', myfontsize,'Fontweight','bold')
yticks([])
%ylabel('cDR')
ax = gca;
ax.XAxis.FontSize = myfontsize;
ax.YAxis.FontSize = myfontsize;
sgtitle('cDR, 3 separate experiments','FontSize', myfontsize,'Fontweight','bold');
%
% too messy, forget it:

% now DR, not cDR

% extract each experiment separately
% the rows are indexed by the concentrations
% the columns are times 1 through 8

exp1 = [exp11' exp12' exp13' exp14' exp15' exp16' exp17' exp18'];
exp2 = [exp21' exp22' exp23' exp24' exp25' exp26' exp27' exp28'];
exp3 = [exp31' exp32' exp33' exp34' exp35' exp36' exp37' exp38'];

% next do TIME derivatives to show adaptation - bad that estimates include negative numbers

%
% derivative could be computed this way:
% using function cent_diff_n(f,h,n) which computes 3-point central
% difference (with n=3) and spep size (h=1) of f. On boundary points
% it uses best estimate with less points
%dexp1 = zeros(noconcentrations,8);
%for i=1:noconcentrations
%    dexp1(i,:) = cent_diff_n(exp1(i,:),1,3);
%end
%% transpose so that now each row is a time and the columns are concentrations
%dexp1t = dexp1';

% but I don't like that the initial value of y is then estimated to nonzero, which is not realistic; so go back to plain diff and add a zero at the start

figure(3)
t = tiledlayout(1, 3, "TileSpacing", "tight");
nexttile
dexp1 = [zeros(noconcentrations,1) diff(exp1')'];
plot(dexp1','linewidth',mylinewidth)
%legend('input 1','input 2','input 3','input 4','input 5','input 6','input 7','input 8','FontSize', myfontsize,'location','northeast','Fontweight','bold')
xlabel('times','Fontweight','bold')
ylabel('estimated output','Fontweight','bold');
ax = gca;
ax.XAxis.FontSize = myfontsize;
ax.YAxis.FontSize = myfontsize;
xticks(1:8)
xlim([1 8])
%yticks([])
%axis square

nexttile
dexp2 = [zeros(noconcentrations,1) diff(exp2')'];
plot(dexp2','linewidth',mylinewidth)
legend('input 1','input 2','input 3','input 4','input 5','input 6','input 7','input 8','FontSize', myfontsize,'location','northeast','Fontweight','bold')
xlabel('times','Fontweight','bold')
%ylabel('estimated output','Fontweight','bold');
ax = gca;
ax.XAxis.FontSize = myfontsize;
ax.YAxis.FontSize = myfontsize;
xticks(1:8)
xlim([1 8])
yticks([])
%axis square

nexttile
dexp3 = [zeros(noconcentrations,1) diff(exp3')'];
plot(dexp3','linewidth',mylinewidth)
%legend('input 1','input 2','input 3','input 4','input 5','input 6','input 7','input 8','FontSize', myfontsize,'location','northeast','Fontweight','bold')
xlabel('times','Fontweight','bold')
%ylabel('estimated output','Fontweight','bold');
ax = gca;
ax.XAxis.FontSize = myfontsize;
ax.YAxis.FontSize = myfontsize;
yticks([])
xticks(1:8)
xlim([1 8])
%axis square
sgtitle('Adaptation','FontSize', myfontsize,'Fontweight','bold');



