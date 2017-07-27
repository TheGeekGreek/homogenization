f = @(x) exp(10 * x);
A_epsilon = @(x,epsilon) 2 + cos(2 * pi * x/epsilon);
u_0 = @(x) 1/sqrt(3) * 1/100 * ((exp(10) - 1 ) * x + 1 - exp(10 * x));
u_0_der = @(x) 1/sqrt(3) * 1/100 * (exp(10) - 10 * exp(10 * x));
psi_epsilon = @(x,epsilon) min(1,min(x,1-x)/epsilon);
N_der = @(x) 1 - sqrt(3)./(2 + cos(2 * pi * x));
N = @(x) integral(N_der,0,x);
w_1_epsilon = @(x,epsilon) u_0(x) - epsilon * psi_epsilon(x,epsilon) ...
      * N(x/epsilon) * u_0_der(x);
x = linspace(0,1);
epsilon = .01;
y_1 = arrayfun(@(x) w_1_epsilon(x,epsilon),x);

c = 1/10 * integral(@(t) f(t)./A_epsilon(t,epsilon),0,1) * ...
    1/integral(@(t) 1./A_epsilon(t,epsilon),0,1);
u_epsilon = @(x,epsilon) integral(@(t) 1./A_epsilon(t,epsilon) * ...
    (c - 1/10 * exp(10 * t)),0,x,'ArrayValued',true);
y_2 = arrayfun(@(x) u_epsilon(x,epsilon),x);

plot(x,y_1)
hold on;
plot(x,y_2)
grid on;
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 14);
l = legend('$u_\varepsilon$','$w_\varepsilon^1$');
set(l, 'fontsize', 14, 'interpreter', 'latex','Location','northwest');
fig = gcf;
fig.PaperUnits = 'centimeters';