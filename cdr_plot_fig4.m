%  for paper, removing not-used stuff from cdr_iffl_plot.m

function  cdr_plot_fig4

myfontsize = 36;
mylinewidth = 4;
step = 0.01;
showplots = 1;
tf1 = 1;
Ulow = 0.01
Umax = 100
nsteps=100
x0 = 0
y0 = 1/10
gain = 10

tspan1 = linspace(0,tf1,nsteps);

howmanystates = 3;

U = Ulow:step:Umax;
I = length(U);
i=0;

allresults=[];
allresultsy=[];
for u = U
  i = i+1;
 fprintf('input value %d out of %d\n',i,I)
  W0 = [x0 y0 0];
  options = odeset('RelTol',1e-10, 'AbsTol',1e-12);
%  [ts,x] = ode45(@f,tspan1,W0,options);
  [ts,x] = ode23s(@f,tspan1,W0,options);
  y=x(end,2); 
  z=x(end,3);
  allresults=[allresults;u z];
  allresultsy=[allresultsy;u y];
end
ALLRESULTS = allresults;
ALLRESULTSY = allresultsy;

if showplots==1
 figure(1)
 t = tiledlayout(1, 3, "TileSpacing", "tight");
 nexttile
 plot(allresultsy(:,1),allresultsy(:,2),'LineWidth',mylinewidth)
 xlabel('input','Fontweight','bold')
 ylabel('DR','Fontweight','bold');
 set(gca, 'XScale', 'log')
 ax = gca;
 ax.XAxis.FontSize = myfontsize;
 ax.YAxis.FontSize = myfontsize;
 axis square
 set(gca,'XTick',[0.1  1E1 1E2 1E3])

 nexttile
 plot(allresults(:,1),allresults(:,2),'LineWidth',mylinewidth)
 xlabel('input','Fontweight','bold')
 ylabel('cDR','Fontweight','bold');
 set(gca, 'XScale', 'log')
 ax = gca;
 ax.XAxis.FontSize = myfontsize;
 ax.YAxis.FontSize = myfontsize;
 axis square
 set(gca,'XTick',[0.1  1E1 1E2 1E3])
 sgtitle('DR and cDR for IFFL1 example','FontSize',myfontsize,'Fontweight','bold');
end

% -----------------------------------------------------------------------
  function dxdt = f(t,x)
      
  dxdt = zeros(howmanystates,1); 

  X=x(1);
  Y=x(2);
  Z=x(3);

  dxdt(1) = -X+u;
  dxdt(2) = - gain*X*Y + u;
  dxdt(3) = Y;
  end
end
