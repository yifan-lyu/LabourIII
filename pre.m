%% Parenting as investment (Cunha 2007), Yifan Lyu, SSE
clear; clc;


syms delta_t gamma1 theta_t phi_t gamma2 I_t phi gamma3 theta_p rho_t

%theta_1 = 1;     % initial draw of human capital
theta_plus = delta_t*(gamma1*theta_t.^phi_t + gamma2*I_t.^phi_t + gamma3*theta_p.^phi_t).^(rho_t/phi_t);
%showlatex(theta_plus);

increment1_d = diff(theta_plus,I_t);
increment2_d = diff(increment1_d,theta_t);

showlatex(increment1_d);




function showlatex(equ_price)
% this function converts symbolic expression to latex expression
% Yifan Lyu, 2020
LaTeX_Expr = latex(equ_price);
display(LaTeX_Expr);
clipboard('copy', LaTeX_Expr);
text(0.5, 0.5, ['$$' LaTeX_Expr '$$'], 'Interpreter','latex', 'FontSize',30, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle')
axis off;
end




