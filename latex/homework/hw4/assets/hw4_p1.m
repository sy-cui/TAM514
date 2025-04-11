clear; close all;

%% 
mu = linspace(0, 20, 128);
func = @(x) 1+cos(x)*cosh(x);
guesses = pi*(0.5+0:6);
roots = guesses;
for i=1:length(guesses)
    roots(i) = fsolve(func, guesses(i));
end

figure
    plot(mu, -cos(mu), '-k', LineWidth=1.5, ...
        DisplayName="$-\cos\mu$")
    hold on
    plot(mu, 1./cosh(mu), '--k', LineWidth=1.5, ...
        DisplayName="$(\cosh\mu)^{-1}$")
    yline(0, '-.k', LineWidth=1, HandleVisibility='off')
    scatter(roots, -cos(roots), 60, 'k', 'filled', HandleVisibility='off')
    hold off
    legend(Interpreter="latex", location="northeast")
    xlabel("$\mu$", Interpreter="latex")
    xlim([0, 20])
    set(gca, Fontsize=20, Fontname="Times new roman")
papersize = [540 360];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw4_p1_evals.pdf -bestfit
