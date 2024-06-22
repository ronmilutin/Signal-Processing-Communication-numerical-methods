a = 1;
b = 5;
s = 3^(1/4);
tol = 10^(-12); 
ID1 = 319095832;
ID2 = 316389584; 
Xo = a + (b-a)*(ID1/(ID1+ID2)); 
format long 

%---Q1---%

[X_n_1,err_n,Xn_Diff_list_1] = Newton_Raphson(Xo,tol,s,'Q1'); 

%---B---%
n = [1:length(err_n)]'; 
X_n = X_n_1(:,1:length(X_n_1)-1)'; 
X_n_Diff = Xn_Diff_list_1'; 
Error_n = err_n'; 
Table1 = table(n,X_n, X_n_Diff, Error_n); 
disp(Table1); 
fprintf('It took %d iteration to converge\n',length(n));

%---c---%
figure(1)
x_axis = log(err_n(1:end-1)); 
y_axis = log(err_n(2:end));  
plot(x_axis,y_axis,'--o');
title('Newton Raphson');
xlabel('log(\bf\epsilon_n_-_1)'); 
ylabel('log(\bf\epsilon_n)');
grid on;
movegui(figure(1),"north")

%---Q2---%

X1 = Xo + (b-Xo)*(ID1/(ID1+ID2)); 
[X_n_2,err_n,Xn_Diff_list_2] = Secant_method(Xo,X1,tol,s); 

%---A---%

n = [1:length(err_n)]'; 
X_n = X_n_2(:,1:length(X_n_2)-1)'; 
X_n_Diff = Xn_Diff_list_2'; 
Error_n = err_n'; 
T2 = table(n,X_n, X_n_Diff, Error_n); 
disp(T2);
fprintf('It took %d iteration to converge\n',length(n));

figure(2)
x_axis = log(err_n(1:end-1)); 
y_axis = log(err_n(2:end)); 
plot(x_axis,y_axis,'--o');
title('Secant Method');
xlabel('log(\bf\epsilon_n_-_1)');
ylabel('log(\bf\epsilon_n)');
grid on;
movegui(figure(2),"south")

%---Q3---%

%---A---%

Xo = 5;
s3 = 2;
[X_n_3a,err_n,Xn_Diff_list_3_A] = Newton_Raphson_multiple_sqrt(Xo,1,tol,s3,'f'); 

n = [1:length(err_n)]'; 
X_n = X_n_3a(:,1:length(X_n_3a)-1)'; 
X_n_Diff = Xn_Diff_list_3_A'; 
Error_n = err_n'; 
T3_A = table(n,X_n, X_n_Diff, Error_n); 
fprintf('%60s\n','<strong> NR0 </strong>');
disp(T3_A); 
fprintf('It took %d iteration to converge\n',length(n));

figure(3)
subplot(3,1,1);
x_axis = log(err_n(1:end-1)); 
y_axis = log(err_n(2:end)); 
plot(x_axis,y_axis,'--o');
title('Newton Raphson Multiple Sqrt');
xlabel('log(\bf\epsilon_n_-_1)');
ylabel('log(\bf\epsilon_n)');
grid on;

%---B---%

[X_n_3_B,err_n,Xn_Diff_list_3_B] = Newton_Raphson_multiple_sqrt(Xo,1,tol,s3,'u'); 

n = [1:length(err_n)]'; 
X_n = X_n_3_B(:,1:length(X_n_3_B)-1)'; 
X_n_Diff = Xn_Diff_list_3_B';
Error_n = err_n'; 
T3_B = table(n,X_n, X_n_Diff, Error_n); 
fprintf('%60s\n','<strong> NR1 </strong>');
disp(T3_B); 
fprintf('It took %d iteration to converge\n',length(n));

subplot(3,1,2);
x_axis = log(err_n(1:end-1)); 
y_axis = log(err_n(2:end)); 
plot(x_axis,y_axis,'--o');
title('Newton Raphson Multiple Sqrt');
xlabel('log(\bf\epsilon_n_-_1)');
ylabel('log(\bf\epsilon_n)');
grid on;


%---C---%

[X_n_3_C,err_n,Xn_Diff_list_3_C] = Newton_Raphson_multiple_sqrt(Xo,2.999,tol,s3,'f'); 

n = [1:length(err_n)]'; 
X_n = X_n_3_C(:,1:length(X_n_3_C)-1)'; 
X_n_Diff = Xn_Diff_list_3_C'; 
Error_n = err_n'; 
T3_C = table(n,X_n, X_n_Diff, Error_n); 
fprintf('%60s\n','<strong> NR2 </strong>'); 
disp(T3_C); 
fprintf('It took %d iteration to converge\n',length(n));

subplot(3,1,3);
x_axis = log(err_n(1:end-1)); 
y_axis = log(err_n(2:end)); 
plot(x_axis,y_axis,'--o');
title('Newton Raphson Multiple Sqrt');
xlabel('log(\bf\epsilon_n_-_1)');
ylabel('log(\bf\epsilon_n)');
grid on;
movegui(figure(3),"west")


%---Q4---%

Xo = pi/2;
s4 = 1.895494267034;
[X_n_4_A,err_n,Xn_Diff_list_4_A] = Fixed_Point(Xo,tol,s4,'A'); 

n = [1:length(err_n)]'; 
X_n = X_n_4_A(:,1:length(X_n_4_A)-1)'; 
X_n_Diff = Xn_Diff_list_4_A';
Error_n = err_n'; 
T4_A = table(n,X_n, X_n_Diff, Error_n); 
disp(T4_A); 
fprintf('It took %d iteration to converge\n',length(n));

figure(4)
subplot(3,1,1);
x_axis = log(err_n(1:end-1)); 
y_axis = log(err_n(2:end));
plot(x_axis,y_axis,'--o');
title('Fixed Point Method');
xlabel('log(\bf\epsilon_n_-_1)');
ylabel('log(\bf\epsilon_n)');
grid on;

%---B---%

[X_n_4_B,err_n,Xn_Diff_list_4_B] = Newton_Raphson(Xo,tol,s4,'Q4'); 

n = [1:length(err_n)]'; 
X_n = X_n_4_B(:,1:length(X_n_4_B)-1)'; 
X_n_Diff = Xn_Diff_list_4_B'; 
Error_n = err_n'; 
T4_B = table(n,X_n, X_n_Diff, Error_n);
disp(T4_B); 
fprintf('It took %d iteration to converge\n',length(n));

subplot(3,1,2);
x_axis = log(err_n(1:end-1)); 
y_axis = log(err_n(2:end));
plot(x_axis,y_axis,'--o');
title('Question 4 Part B: Newton Raphson Method');
xlabel('log(\bf\epsilon_n_-_1)');
ylabel('log(\bf\epsilon_n)');
grid on;


%---D---%

s4d=0;
Xo = 1/2;
[X_n_4_D,err_n,Xn_Diff_list_4_D] = Fixed_Point(Xo,tol,s4d,'D'); 

n = [1:length(err_n)]'; 
X_n = X_n_4_D(:,1:length(X_n_4_D)-1)'; 
X_n_Diff = Xn_Diff_list_4_D'; 
Error_n = err_n'; 
T4_D = table(n,X_n, X_n_Diff, Error_n); 
disp(T4_D);
fprintf('It took %d iteration to converge',length(n));

subplot(3,1,3);
x_axis = log(err_n(1:end-1)); 
y_axis = log(err_n(2:end)); 
plot(x_axis,y_axis,'--o');
title('Fixed Point Method_D');
xlabel('log(\bf\epsilon_n_-_1)');
ylabel('log(\bf\epsilon_n)');
grid on;
movegui(figure(4),"east");

%---Func---%

%---Newton_Raphson---%

function [X_n,err_n,Xn_Diff_list] = Newton_Raphson(Xo,tol,s,func)
    iter = 2; 
    err_n = zeros;  
    Xn_Diff_list = zeros; 
    X_n(1) = Xo; 
    X_n(2) = Xo - h(Xo,func);
    Xn_Diff = abs(X_n(iter)-X_n(iter-1)); %|x_(n)-x_(n-1)|
    Xn_Diff_list(1) = Xn_Diff;
    err_n(1) = abs(X_n(1) - s); %|x_n - s|
    while Xn_Diff >= tol                     
        X_n(iter+1) = X_n(iter) - h(X_n(iter),func);
        err_n(iter) = abs(X_n(iter) - s);     
        iter=iter+1;    
        Xn_Diff = abs(X_n(iter)-X_n(iter-1));  
        Xn_Diff_list(iter-1) = Xn_Diff;       
    end
end

function fx = my_func(x,func) %Calc f(x) depended on function needed.
    if func == 'A'
        fx = x^4-3;
    elseif func == 'B'
        fx = x-2*sin(x);
    end
end

function f_tag_x = derivative_func(x,func) %Calc f'(x) depended on function needed.
    if func == 'A'
        f_tag_x = 4*x^3;
    elseif func == 'B'
        f_tag_x = 1-2*cos(x);
    end
end

function div = h(x,func)  %Calc h(x) = f(x) / f'(x)
    if func == 'Q1'
        div = my_func(x,'A') / derivative_func(x,'A');
    elseif func == 'Q4'
        div = my_func(x,'B') / derivative_func(x,'B'); 
    end
        
end


%---Newton_Raphson_multiple_sqrt---%


function [X_n,err_n,Xn_Diff_list] = Newton_Raphson_multiple_sqrt(Xo,mult,tol,s,func)
    iter = 2; 
    err_n = zeros; 
    Xn_Diff_list = zeros;
    X_n(1) = Xo;
    X_n(2) = Xo - h_2(Xo,func);
    Xn_Diff = abs(X_n(iter)-X_n(iter-1));
    Xn_Diff_list(1) = Xn_Diff;
    err_n(1) = abs(X_n(1) - s);
    while (abs(X_n(iter)-X_n(iter-1))>tol) && (abs(X_n(iter)-s)<(abs(X_n(iter-1)-s))) 
        X_n(iter+1) = X_n(iter) - mult * h_2(X_n(iter),func);
        err_n(iter) = abs(X_n(iter) - s);   
        iter=iter+1;    
        Xn_Diff = abs(X_n(iter)-X_n(iter-1));  
        Xn_Diff_list(iter-1) = Xn_Diff; 
    end
    
end

function fx = my_func_2(x)  %f(x)
    fx = x^5-6*x^4+14*x^3-20*x^2+24*x-16;
end

function f_tag_x = derivative_func_2(x) %f'(x)
    f_tag_x = 5*x^4-24*x^3+42*x^2-40*x+24;
end

function f_tag_tag_x = deri_deri_func_2(x) %f''(x)
    f_tag_tag_x = 20*x^3-72*x^2+84*x-40;
end

function y = u(x) % u = f(x) / f'(x) ;
    y = my_func_2(x) / derivative_func_2(x);
end

function y = deri_u(x) % u'(x) = 1 - (f(x) * f''(x))/(f'(x) * f'(x))
    y =  1-(my_func_2(x)*deri_deri_func_2(x))/(derivative_func_2(x)^2);
end

function y = h_u(x) % h = u(x) / u'(x)
        y = u(x)/deri_u(x);
end

function [y] = h_2(x,func) 
    if func == 'f'
        y = my_func_2(x) / derivative_func_2(x);
    elseif func == 'u'
        y = h_u(x);
    end
    
end

%---Fixed_Point---%

function [X_n,err_n,Xn_Diff_list] = Fixed_Point(x0,tol,s,func)
    iter = 2;
    err_n = zeros;
    Xn_Diff_list = zeros;
    X_n(1) = x0;
    X_n(2) = g(X_n(1), func);
    while abs(X_n(iter) - X_n(iter-1)) >= tol
        X_n(iter+1) = g(X_n(iter), func);
        err_n(iter) = abs(X_n(iter) - s);     
        iter=iter+1;    
        Xn_Diff = abs(X_n(iter)-X_n(iter-1));  
        Xn_Diff_list(iter-1) = Xn_Diff;
    end
end

function g = g(x,func) 
    if func == 'A' 
        g = 2*sin(x);
    elseif func == 'D'
        g = asin(x/2);
    end
end

%---Secant_method---%

function [X_n,err_n,Xn_Diff_list] = Secant_method(Xo, x1, tol,s)
    iter = 2; 
    err_n = zeros; 
    Xn_Diff_list = zeros; 
    X_n(1) = Xo;
    X_n(2) = x1;
    Xn_Diff_list(1) = abs(X_n(2)-X_n(1));
    err_n(1) = abs(X_n(1) - s); 
    while abs(X_n(iter) - X_n(iter-1)) > tol
        X_n(iter+1) = X_n(iter) - func(X_n(iter))*(X_n(iter)-X_n(iter-1))/(func(X_n(iter))-func(X_n(iter-1))); 
        err_n(iter) = abs(X_n(iter) - s);  
        iter=iter+1;    
        Xn_Diff_list(iter-1) = abs(X_n(iter)-X_n(iter-1));
    end
end

function fx = func(x) 
    fx = x^4-3; 
end