% for paper; was example1.m; inserting cdr_integral

function cdr_plot_fig5

global howmanystates p u K

myfontsize = 36;
mylinewidth = 4;

ulow = 10;
uhigh = 65;

x0 = 0.1
p = 6
K = 0

figure(1)

display ('computing first time')

[PARAMS3,VALUES3,ALLRESULTS3] = cdr_integral([ulow  uhigh    x0    p p     K    3]);VALUES3
plot(ALLRESULTS3(:,1),ALLRESULTS3(:,2),'LineWidth',mylinewidth)
hold on

display ('computing second time')

[PARAMS4,VALUES4,ALLRESULTS4] = cdr_integral([ulow  uhigh    x0    p p     K    4]);VALUES4
plot(ALLRESULTS4(:,1),ALLRESULTS4(:,2),'LineWidth',mylinewidth)

display ('computing third time')

[PARAMS5,VALUES5,ALLRESULTS5] = cdr_integral([ulow  uhigh    x0    p p     K    5]);VALUES5
plot(ALLRESULTS5(:,1),ALLRESULTS5(:,2),'LineWidth',mylinewidth)

display ('computing fourth time')

[PARAMS6,VALUES6,ALLRESULTS6] = cdr_integral([ulow  uhigh    x0    p p     K    6]);VALUES6
plot(ALLRESULTS6(:,1),ALLRESULTS6(:,2),'LineWidth',mylinewidth)

hold off
set(gca, 'XScale', 'log')
legend('T=3','T=4', 'T=5', 'T=6','FontSize', myfontsize,'location','southeast','Fontweight','bold')
xlabel('input','Fontweight','bold')
ylabel('cDR','Fontweight','bold');
title('cDR','FontSize', myfontsize,'Fontweight','bold');
ax = gca;
ax.XAxis.FontSize = myfontsize;
ax.YAxis.FontSize = myfontsize;
axis square

% -------------------------------------------

function  [PARAMS,VALUES,ALLRESULTS] = cdr_integral(PARAMS)

uselog = 1;
step = 0.1;
showplots = 0;
%PARAMS = [Ulow, Umax, x0,y0,p,K,tf1];

Ulow  = PARAMS(1);
Umax  = PARAMS(2);
x0    = PARAMS(3);
y0    = PARAMS(4);
p     = PARAMS(5);
K     = PARAMS(6);
tf1   = PARAMS(7);
  
% how many intermediate steps to store (more is better for plotting):
nsteps=100;

tspan1 = linspace(0,tf1,nsteps);

howmanystates = 3;

U = Ulow:step:Umax;
I = length(U);
i=0;

if uselog == 1
allresults=[];
allresultsy=[];
for u = U
  i = i+1;
  fprintf('input value %d out of %d\n',i,I)
  W0 = [log(u/(K+x0)) y0-p 0];
  [ts,x] = ode45(@flog,tspan1,W0);% , options);
  y=x(end,2)+p; 
  z=x(end,3);
  if ts(end) < tf1
    errorflag = 1 % ode solution ended before tf
  else
    errorflag = 0;
  end
  if errorflag == 0
    allresults=[allresults;u z];
    allresultsy=[allresults;u y];
  end
end
end % uselog

ALLRESULTS = allresults;
V = allresults(:,2)';
D = diff(V);
S = sign(D);
Z = sign(diff(S));
ZeroCrossings = find(Z ~= 0);
VALUES = [V(1),V(ZeroCrossings),V(end)];

end

% -----------------------------------------------------------------------
 function dxdt = flog(t,x) % using transformation w = log(u/x)

  dxdt = zeros(howmanystates,1);    % a column vector

  X=x(1);
  Y=x(2);
  Z=x(3);

  % initial condition is log(u/(K+x0)),0,0; integrating original 'y', which is Y+p
  dxdt(1) = - (1-K*exp(X)/u)*Y;
  dxdt(2) = exp(X) - (Y+p);
  dxdt(3) = Y+p;

 end

% -----------------------------------------------------------------------

end
