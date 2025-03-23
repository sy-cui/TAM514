clear; close all;

%% 
t0 = linspace(-pi/2+0.1, pi/2-0.1, 128);
y = tan(t0);
t1 = t0 + pi;
t2 = t1 + pi;
t3 = t2 + pi;
x1 = [pi, 2*pi, 3*pi];
x2 = [0.5*pi, 1.5*pi, 2.5*pi];

k = 0.8;
tx = linspace(0, 10, 2);
eta1 = fsolve(@(x) tan(x) + k*x, pi);
eta2 = fsolve(@(x) tan(x) + k*x, 2*pi);
eta3 = fsolve(@(x) tan(x) + k*x, 3*pi);
etas = [eta1, eta2, eta3];

figure
hold on 
xline([pi/2, 3*pi/2, 5*pi/2], '--k', HandleVisibility='off')
yline(0, '--k', HandleVisibility='off')

plot(t0, y, '-k', Linewidth=1.5, HandleVisibility='off')
plot(t1, y, '-k', Linewidth=1.5, HandleVisibility='off')
plot(t2, y, '-k', Linewidth=1.5, HandleVisibility='off')
plot(t3, y, '-k', Linewidth=1.5, HandleVisibility='off')
plot(tx, -k*tx, '-.k', Linewidth=1.5, HandleVisibility='off')

scatter(x1, 0*x1, 100, '^k', MarkerFaceColor='k', MarkerFaceAlpha=0.2, ...
    DisplayName="Fixed-fixed $\eta_j$")
scatter(x2, 0*x2, 100, 'ok', MarkerFaceColor='k', MarkerFaceAlpha=0.2, ...
    DisplayName="Fixed-free $\eta_j$")
scatter(etas, tan(etas), 50, 'square', MarkerEdgeColor='k', ...
    MarkerFaceColor='k',  DisplayName="$\eta_j$")

hold off; box on;
xlim([0, 3*pi+0.5]);
ylim([-10, 10]);
xlabel("$\eta$", Interpreter='latex')
legend(Interpreter='latex', Location="northoutside", NumColumns=3)
set(gca, Fontsize=20, Fontname="Times New Roman");
papersize = [540 300];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw1_p3_eval.eps -bestfit