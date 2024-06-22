clear all
close all
clc

%---Q1---
%---A----

theta_sample = [2,3,4,5]; 
theta_restore = linspace(0,pi,41); 
for j = 1:length(theta_sample)
    arbitrary_theta = linspace(0,pi,theta_sample(j)); 
    for i = 1:length(theta_restore)
        Phi_approx(j,i) = Lagrange_Interpolation(theta_restore(i),arbitrary_theta,'A');
    end
end
for i = 1:length(theta_restore)  
   Phi_real(i) = potential(theta_restore(i),'A');
end
figure(1)
subplot(3,2,1)
p = plot(theta_restore, Phi_real, theta_restore, Phi_approx,'--.');
p(1).LineWidth = 1.3;
legend('Real Phi','Approx Phi n+1=2','Approx Phi n+1=3','Approx Phi n+1=4','Approx Phi n+1=5','Location','northwest')
title('Q1-A- \phi(\theta)')
xlabel("\theta")
ylabel("\phi(\theta)")
grid on
movegui('northeast')

%---B---

relative_err = 2:2:20;
for i = 1:length(relative_err)
    rel_err_B(i) = Relative_Error(theta_restore,relative_err(i),'A');
end
subplot(3,2,2)
semilogy(relative_err,rel_err_B,'--o',LineWidth=1.3)
title('Q1-B- Relative Error')
xlabel("n+1")
ylabel("Relative Error")
grid on


%---C---

theta_sample = [3,7,11,15]; 
theta_restore = linspace(0,pi,41);
for j = 1:length(theta_sample)
    arbitrary_theta = linspace(0,pi,theta_sample(j)); 
    for i = 1:length(theta_restore)
        Phi_approx_C(j,i) = Lagrange_Interpolation(theta_restore(i),arbitrary_theta,'C');
    end
end
for i = 1:length(theta_restore)  
   Phi_real_C(i) = potential(theta_restore(i),'C');
end
subplot(3,2,3)
p = plot(theta_restore, Phi_real_C, theta_restore, Phi_approx_C,'--.',LineWidth=1.3);
legend('Real Phi','Approx Phi:n+1=3','Approx Phi:n+1=7','Approx Phi:n+1=11','Approx Phi:n+1=15','Location','northeast')
title('Q1-C - \phi(\theta)')
xlabel("\theta")
ylabel("\phi(\theta)")
grid on


%----rel_error---
for i = 1:length(relative_err)
    rel_err_C(i) = Relative_Error(theta_restore,relative_err(i),'C');
end
subplot(3,2,4)
semilogy(relative_err,rel_err_C,'--o',LineWidth=1.3)
title('Q1-C - RelativeError')
xlabel("n+1")
ylabel("Relative Error")
grid on


%---D---

theta_sample = [3,7,11,15]; 
theta_restore = linspace(0,pi,41); 
for j = 1:length(theta_sample)
    for i = 1:length(theta_restore)
        Phi_approx_D(j,i) = Lagrange_Interpolation(theta_restore(i),Chebyshev_Roots(theta_sample(j)-1,0,pi),'C');
    end
end
for i = 1:length(theta_restore)  
   Phi_real_D(i) = potential(theta_restore(i),'C');
end
subplot(3,2,5)
p = plot(theta_restore, Phi_real_D, theta_restore, Phi_approx_D,'--.',LineWidth=1.3);
legend('Real Phi','Approx Phi:n+1=3','Approx Phi:n+1=7','Approx Phi:n+1=11','Approx Phi:n+1=15','Location','northeast')
title('Q1-D- \phi(\theta) Chebyshev')
xlabel("\theta")
ylabel("\phi(\theta)")
grid on


%---rel_error---

for i = 1:length(relative_err)
    rel_err_D(i) = Relative_Error_2(theta_restore,relative_err(i),'C');
end
subplot(3,2,6)
semilogy(relative_err,rel_err_D,'--o',relative_err,rel_err_C,'--o',LineWidth=1.3)
legend('Relative error_with chebyshev','Relative error_Uniform Distribution')
title('Q1-D- RelativeError')
xlabel("n+1")
ylabel("Relative Error")
grid on



%---Q2---
%----B----

theta_sample = [2,3,4];
for j = 1:length(theta_sample)
    arbitrary_theta = linspace(0,pi,theta_sample(j));
    for i = 1:length(arbitrary_theta)
        y_1(i) = potential(arbitrary_theta(i),'2');
    end
    [a(j),b(j),c(j)] = find_abc(arbitrary_theta,y_1);
    Phi_LS_2B(j,:) = a(j)+b(j)*sin(theta_restore)+c(j)*cos(theta_restore);
end
for i = 1:length(theta_restore)
    Phi_real_2B(i) = potential(theta_restore(i),'2');
end
figure(2)
subplot(2,2,1)
p = plot(theta_restore, Phi_real_2B, theta_restore, Phi_LS_2B,'--.',LineWidth=1.3);
title("Q2-B- \phi(\theta) Least Squares")
xlabel("\theta")
ylabel("\phi(\theta)")
legend("Real Phi", "Approx Phi:n+1=2","Approx Phi:n+1=3","Approx Phi:n+1=4")
grid on
movegui('north')

%---C---

r_0 = 10;
r=r_0*2.^(0:-1:-8);
theta_sample = 4;
arbitrary_theta = linspace(0,pi,theta_sample);
for j = 1:length(r)
    for i = 1:length(arbitrary_theta)
        y_1(i) = potential_2(arbitrary_theta(i),r(j));
    end
    [a(j),b(j),c(j)] = find_abc(arbitrary_theta, y_1);
    Phi_LS_2C(j,:) = a(j)+b(j)*sin(theta_restore)+c(j)*cos(theta_restore);
    for i = 1:length(theta_restore)
        Phi_real_2C(j,i) = potential_2(theta_restore(i),r(j));
    end
end
relative_error_2_C = sqrt(sum((Phi_LS_2C-Phi_real_2C).^2,2))./sqrt(sum(Phi_real_2C.^2,2));
subplot(2,2,2)
loglog(r', relative_error_2_C, "--o",LineWidth=1.3)
title("Q2-C- RelativeError")
xlabel("Radius[m]")
ylabel("relative error")
grid on


%---D---

theta_sample = 2.^(2:18);
for j = 1:length(theta_sample)
    arbitrary_theta = linspace(0,pi,theta_sample(j)); 
    for i = 1:length(arbitrary_theta)
        y(i) = potential(arbitrary_theta(i),'2');
    end
    y_error = (1+(rand(1,theta_sample(j))-0.5)*10^-1).*y;
    [a(j),b(j),c(j)] = find_abc(arbitrary_theta, y_error);
    Phi_LS_2_D(j,:) = a(j)+b(j)*sin(theta_restore)+c(j)*cos(theta_restore);
end
for i = 1:length(theta_restore)
    Phi_real_2_D(i) = potential(theta_restore(i),'2');
end
relative_error_2_D = sqrt(sum((Phi_LS_2_D-Phi_real_2_D).^2,2))./sqrt(sum(Phi_real_2_D.^2,2));
subplot(2,2,3)
loglog(theta_sample', relative_error_2_D, '--o',LineWidth=1.3)
title("Q2-D- RelativeError")
xlabel("n+1")
ylabel("relative error")
grid on


%---D2---

theta_sample = 2.^(2:18);
for j = 1:length(theta_sample)
    arbitrary_theta = linspace(0,pi,theta_sample(j)); 
    for i = 1:length(arbitrary_theta)
        y_1(i) = potential(arbitrary_theta(i),'2');
    end
    y_error_2 = (1+(rand(1,theta_sample(j))-0.5)*10^-4).*y_1;
    [a(j),b(j),c(j)] = find_abc(arbitrary_theta, y_error_2);
    Phi_LS_2_D_2(j,:) = a(j)+b(j)*sin(theta_restore)+c(j)*cos(theta_restore);
end
for i = 1:length(theta_restore)
    Phi_real_2_D_2(i) = potential(theta_restore(i),'2');
end
relative_error_2_D_2 = sqrt(sum((Phi_LS_2_D_2-Phi_real_2_D_2).^2,2))./sqrt(sum(Phi_real_2_D_2.^2,2));
subplot(2,2,4)
loglog(theta_sample', relative_error_2_D_2, '--o',LineWidth=1.3)
title("Q2-D- Relative Error-Second Delta")
xlabel("n+1")
ylabel("relative error")
grid on


%----Q3----
%----A-----

clear
format long
a = 0;
b = 1;
Integral_Trapezoid = Trapezoid_Integration(a,b);
Integral_Simpson = Simpson_Integration(a,b);
Integral_Real = 4/pi*atan(1);
Err_Trapezoid = abs(Integral_Real - Integral_Trapezoid);
Err_Simpson = abs(Integral_Real - Integral_Simpson);
disp('Trapezoid Error value: ')
disp(Err_Trapezoid)
disp('Trapezoid Integral value: ')
disp(Integral_Trapezoid)
disp('Simpson Error value: ')
disp(Err_Simpson)
disp('Simpson Integral value: ')
disp(Integral_Simpson)

%----B-----

n_list = [5 9 17 33 65 129 257 513];
a=0;
b=pi;
Integral_Simp1 = [];
Integral_Simp2 = [];
Integral_Simp3 = []; 
Integral_Trap1 = [];
Integral_Trap2 = [];
Integral_Trap3 = [];
err_Simp1 = [];
err_Simp2 = [];
err_Simp3 = [];
err_Trap1 = [];
err_Trap2 = [];
err_Trap3 = [];
for theta_sample = n_list
    [s_a, s_b, s_c] = Simpson_composite_Integration(0, pi, theta_sample);
    [t_a, t_b, t_c] = Trapezoid_composite_Integration(0, pi, theta_sample);
    Integral_Simp1 = [Integral_Simp1 s_a];
    Integral_Simp2 = [Integral_Simp2 s_b];
    Integral_Simp3 = [Integral_Simp3 s_c]; 
    Integral_Trap1 = [Integral_Trap1 t_a];
    Integral_Trap2 = [Integral_Trap2 t_b];
    Integral_Trap3 = [Integral_Trap3 t_c];
end

for i = Integral_Simp1
    err_Simp1 = [err_Simp1 abs((i-Integral_Simp1(end))/Integral_Simp1(end))];
end
for i = Integral_Simp2
    err_Simp2 = [err_Simp2 abs((i-Integral_Simp2(end))/Integral_Simp2(end))];
end
for i = Integral_Simp3
    err_Simp3 = [err_Simp3 abs((i-Integral_Simp3(end))/Integral_Simp3(end))];
end
for i = Integral_Trap1
    err_Trap1 = [err_Trap1 abs(i-Integral_Trap1(end))/abs(Integral_Trap1(end))];
end
for i = Integral_Trap2
    err_Trap2 = [err_Trap2 abs(i-Integral_Trap2(end))/abs(Integral_Trap2(end))];
end
for i = Integral_Trap3
    err_Trap3 = [err_Trap3 abs(i-Integral_Trap3(end))/abs(Integral_Trap3(end))];
end


figure(3)
semilogy(n_list,err_Trap1, 'r--o', n_list,err_Trap2,'b--o',n_list,err_Trap3, 'k--o', n_list,err_Simp1,'r--x',n_list,err_Simp2, 'b--x', n_list,err_Simp3,'k--x',LineWidth=1)
title('Q3-B- Relative Error function of n+1');
xlabel('n+1 Values');
ylabel('Relative Error');
legend('f1=1 - Trapez', 'f2=sin(x) - Trapez', 'f3=cos(x) - Trapez', 'f1=1 - Simpson' , 'f2=sin(x) - Simpson', 'f3=cos(x) - Simpson','position',[0.68 0.15 0.2 0.2]);
grid on;
movegui('northwest')

%---------------------------Functions--------------------------------------

function [L_N] = Lagrange_Interpolation(x,arbitrary_theta,section)
    sum = 0;
    for i = 1:length(arbitrary_theta)
        numerator = 1;
        denominator = 1;
        for j = 1:length(arbitrary_theta)
            if (j~=i)
                numerator = numerator * (x-arbitrary_theta(j)); 
                denominator = denominator * (arbitrary_theta(i)-arbitrary_theta(j)); 
            end
        end
        sum = sum + potential(arbitrary_theta(i),section) * (numerator/denominator); 
    end
    L_N = sum;
end


function [roots] = Chebyshev_Roots(n, a, b)
    n = n + 1;
    chebyshev_theta = [];
    for i = 1:n
        x = cos(pi*(2*i-1)/(2*n));
        t = ((b-a)*x+b+a)/2;
        chebyshev_theta = [chebyshev_theta t];
    end
    roots = chebyshev_theta;
end


function [value] = radius(theta,sign,r)  
    delta = 5 * 10^(-3) ; 
    if sign == '+'
        value = sqrt((r*cos(theta)).^2+(r*sin(theta) - delta/2).^2);   
    elseif sign == '-'
        value = sqrt((r*cos(theta)).^2+(r*sin(theta) + delta/2).^2);
    end
end    

function [phi] = potential(theta,sect) 
    if sect == 'A'
        r = 0.05;
    elseif sect == 'C'
        r = 4*10^(-3);
    elseif sect == '2'
        r = 0.1;
    end
    q_pos = sum([1 9 0 9 5 8 3 2 1 6 3 8 9 5 8 4]);  % ids-319095832,316389584
    q_neg = -sum([3 1 9 0 9 5 8 3 3 1 6 3 8 9 5 8]);  % ids-319095832,316389584
    phi = (q_pos/(4*pi*radius(theta,'+',r))) + (q_neg/(4*pi*radius(theta,'-',r)));
end

function [phi] = potential_2(theta,r) 
    q_pos = sum([1 9 0 9 5 8 3 2 1 6 3 8 9 5 8 4]);  % ids-319095832,316389584
    q_neg = -sum([3 1 9 0 9 5 8 3 3 1 6 3 8 9 5 8]);  % ids-319095832,316389584
    phi = (q_pos/(4*pi*radius(theta,'+',r))) + (q_neg/(4*pi*radius(theta,'-',r)));
end

function [rel_err] = Relative_Error(theta_restore,n,sect)
    arbitrary_theta = linspace(0,pi,n);
    numerator_sum = 0;
    denominator_sum = 0;
    for i = 1:length(theta_restore)
        numerator_sum = numerator_sum + (Lagrange_Interpolation(theta_restore(i),arbitrary_theta,sect) - potential(theta_restore(i),sect))^2;
        denominator_sum = denominator_sum + (potential(theta_restore(i),sect))^2;
    end
    rel_err = sqrt(numerator_sum / denominator_sum);
end

function [rel_err] = Relative_Error_2(theta_restore,n,sect)
    numerator_sum = 0;
    denominator_sum = 0;
    for i = 1:length(theta_restore)
        numerator_sum = numerator_sum + (Lagrange_Interpolation(theta_restore(i),Chebyshev_Roots(n-1,0,pi),sect) - potential(theta_restore(i),sect))^2;
        denominator_sum = denominator_sum + (potential(theta_restore(i),sect))^2;
    end
    rel_err = sqrt(numerator_sum / denominator_sum);
end


function [a,b,c] = find_abc(theta,y)
    f0 = ones(length(theta),1);
    f1 = sin(theta');
    f2 = cos(theta');
    F = [f0 f1 f2];
    vector = (inv(F'*F))*F'*y';
    a = vector(1);
    b = vector(2);
    c = vector(3);
end

function [g_val] = g_x(x)
    g_val = 4 / (pi*(1+x^2));
end

function I = Trapezoid_Integration(a, b)
    h = b - a;
    x_1 = a;
    x_2 = b;
    I = (g_x(x_1)+g_x(x_2)) * (h/2);
end

function I = Simpson_Integration(a, b)
    h = b-a;
    x_1 = a;
    x_2 = (a+b)/2;
    x_3 = b;
    I = (h/6) * (g_x(x_1)+4*g_x(x_2)+g_x(x_3));
end

function [I1, I2, I3] = Simpson_composite_Integration(a, b, n)
    h = (b-a)/(n-1);
    function y = q1(x)
        y = potential(x,'2');
    end
    function y = q2(x)
        y = potential(x,'2')*sin(x);
    end
    function y = q3(x)
        y = potential(x,'2')*cos(x);
    end
    y1 = [];
    y2 = [];
    y3 = [];
    for x = a:h:b
        y1 = [y1 q1(x)];
        y2 = [y2 q2(x)];
        y3 = [y3 q3(x)];
    end
    
    function summ = f(y)
        yN = [];
        y2N = [];
        i = 1;
        while i <= length(y)
            if mod(i, 2) == 0 
                y2N = [y2N y(i)];
            else
                yN = [yN y(i)];
            end
            i = i+1;
        end
        summ = 2*sum(yN)+4*sum(y2N)-y(1)-y(end);
    end
    I1 = h*f(y1)/3;
    I2 = h*f(y2)/3;
    I3 = h*f(y3)/3;
end

function [I1, I2, I3] = Trapezoid_composite_Integration(a, b, n)
    h = (b-a)/(n-1);
    function y = q1(x)
        y = potential(x,'2');
    end
    function y = q2(x)
        y = potential(x,'2')*sin(x);
    end
    function y = q3(x)
        y = potential(x,'2')*cos(x);
    end
    y1 = [];
    y2 = [];
    y3 = [];
    for x = a:h:b
        y1 = [y1 q1(x)];
        y2 = [y2 q2(x)];
        y3 = [y3 q3(x)];
    end
    I1 = (2*sum(y1)-(q1(a)+ q1(b))/2)*h;
    I2 = (2*sum(y2)-(q2(a)+ q2(b))/2)*h;
    I3 = (2*sum(y3)-(q3(a)+ q3(b))/2)*h;
end