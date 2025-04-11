clear; close all;

%% Parameters
k = 100;
E = 1;
I = 1;
m = 1;
M = 1;
L = 1;
b = (E*I/m)^0.25;

%% 
w = linspace(0, 50, 8192);
dw = w*0;
for i=1:length(w)
    mu = sqrt(w(i))*L / (3*b);
    ka = k*b^3/(E*I*w(i)^1.5);

    A = [sin(2*mu) sinh(2*mu) -sin(mu) -sinh(mu) 0 0;
         cos(2*mu) cosh(2*mu) cos(mu) cosh(mu) 0 0;
         -sin(2*mu) sinh(2*mu) sin(mu) -sinh(mu) 0 0;
         cos(2*mu)+ka*sin(2*mu) -cosh(2*mu)+ka*sinh(2*mu) cos(mu) -cosh(mu) -ka 0;
         -sin(2*mu) -sinh(2*mu) 0 0 2-w(i)^2*M/k -1;
         0 0 0 0 -1 1-w(i)^2*M/k];
    
    dw(i) = det(A);
end

%% Plotting
figure
    semilogy(w, abs(dw), '-k', LineWidth=1.5)
    xlabel("$\omega$", Interpreter="latex")
    ylabel("$|\det(A)|$", Interpreter="latex")
    set(gca, Fontsize=20, Fontname="Times new roman")
    grid on
papersize = [540 360];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw4_p3_detA.pdf -bestfit
