%Definitions
f = @(x) exp(10 * x);
A_epsilon = @(x,epsilon) 2 + cos(2 * pi * x/epsilon);
dN = @(x) 1 - sqrt(3)./(2 + cos(2 * pi * x));
N = @(x) integral(dN,0,x);
u_0 = @(x) 1/sqrt(3) * 1/100 * ((exp(10) - 1 ) * x + 1 - f(x));
du_0 = @(x) 1/sqrt(3) * 1/100 * (exp(10) - 1 - 10 * f(x));
ddu_0 = @(x) -1/sqrt(3) * f(x);
psi_epsilon = @(x,epsilon) min(1,min(x,1-x)/epsilon);
dpsi_epsilon = @(x,epsilon) -1/(2 * epsilon) * sign(2 * x - 1) * (1 + ...
    sign(1 - 1/(2 * epsilon) * (1 - abs(2 * x - 1))));
w_1_epsilon = @(x,epsilon) u_0(x) - epsilon * psi_epsilon(x,epsilon) ...
      * N(x/epsilon) * du_0(x);
t_1 = @(x,epsilon) dpsi_epsilon(x,epsilon) * N(x/epsilon) * du_0(x);
t_2 = @(x,epsilon) 1/epsilon * psi_epsilon(x,epsilon) * dN(x/epsilon) * du_0(x);
t_3 = @(x,epsilon) psi_epsilon(x,epsilon) * N(x/epsilon) * ddu_0(x);
dw_1_epsilon = @(x,epsilon) du_0(x) - epsilon * (t_1(x,epsilon) + ...
    t_2(x,epsilon) + t_3(x,epsilon));
c = @(epsilon) integral(@(t) f(t)./A_epsilon(t,epsilon),0,1) * ...
    1/integral(@(t) 1./A_epsilon(t,epsilon),0,1);
u_epsilon = @(x,epsilon,c) 1/10 * integral(@(t) 1./A_epsilon(t,epsilon) * ...
    (c - f(t)),0,x,'ArrayValued',true);
du_epsilon = @(x,epsilon,c) 1/10 * 1/A_epsilon(x,epsilon) * (c - f(x));  

%Parameters
x = linspace(0,1,1e+3);
epsilon = 1e-2;

%Windows
f1 = figure;
f2 = figure;
f3 = figure;

%Plot of the approximation of u_eps by w_1_eps
c = c(epsilon);
y_1 = arrayfun(@(x) u_epsilon(x,epsilon,c),x);
y_2 = arrayfun(@(x) w_1_epsilon(x,epsilon),x);
figure(f1);
plot(x,y_1)
hold on;
plot(x,y_2)
grid on;
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 24);
l1 = legend('$u_\varepsilon$','$w_\varepsilon^1$');
set(l1, 'fontsize', 24, 'interpreter', 'latex','Location','northwest');
fig1 = gcf;
set(fig1,'Units','centimeters');
pos1 = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos1(3), pos1(4)])
print(fig1,'u_eps_w_1_eps_4','-dpdf','-r0')

%Plot of the approximation of du_eps by dw_1_eps
dy_1 = arrayfun(@(x) du_epsilon(x,epsilon,c),x);
dy_2 = arrayfun(@(x) dw_1_epsilon(x,epsilon),x);
figure(f2);
plot(x,dy_1);
hold on;
plot(x,dy_2);
grid on;
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 24);
l2 = legend('$\partial_x u_\varepsilon$','$\partial_x w_\varepsilon^1$');
set(l2, 'fontsize', 24, 'interpreter', 'latex','Location','southwest');
fig2 = gcf;
set(fig2,'Units','centimeters');
pos2 = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos2(3), pos2(4)])
print(fig2,'du_eps_dw_1_eps_4','-dpdf','-r0')

%Error
epsilon_array = linspace(1e-1,8e-3,6);
n = length(epsilon_array);
error = zeros(n);
for i = 1:n
norm_1 = integral(@(x) (u_epsilon(x,epsilon_array(i),c) - w_1_epsilon(x,epsilon_array(i)))^2,0,1,'ArrayValued',true);
norm_2 = integral(@(x) (du_epsilon(x,epsilon_array(i),c) - dw_1_epsilon(x,epsilon_array(i)))^2,0,1,'ArrayValued',true);
norm = sqrt(norm_1 + norm_2);
error(i) = norm;
end

figure(f3);
loglog(epsilon_array,error,'-x');
hold on;
loglog(epsilon_array,error(1)/epsilon_array(1)^.5 * epsilon_array.^.5,'color','red')
grid on;
xlabel('$\varepsilon$', 'interpreter', 'latex', 'fontsize', 24);
l3 = legend('$\|u_\varepsilon - w_\varepsilon^1\|_{H^1(\Omega)}$','$\mathcal{O}(\sqrt{\varepsilon})$');
set(l3, 'fontsize', 24, 'interpreter', 'latex','Location','northwest');
fig3 = gcf;
set(fig3,'Units','centimeters');
pos3 = get(fig3,'Position');
set(fig3,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos3(3), pos3(4)])
print(fig3,'error','-dpdf','-r0')
