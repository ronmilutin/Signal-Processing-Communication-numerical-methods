%---Q1---%

M = 18; 
p = 1 ; 
h_val = [1,2,5,10,20,50];
h_list = h_val.*((pi * p) ./ M);
Rel_Error_b = zeros(1,6)';                            % Rel Error b 
Rel_Error_c = zeros(1,6)';                            % Rel Error c
Rel_Error_d = zeros(1,6)';                            % Rel Error d
K_A = zeros(1,6)';                                    % Condition number vector
q = [3,1,9,0,9,5,8,3,2,3,1,6,3,8,9,5,8,4]';           % vector of id's
for k = 1:6
    
    %---A---%
    h = h_list(k); 
    A = build_matrix(h,M);    
    v = A * q;                                        % potational vector v
    [L,U,P] = lu(A);
    K_A(k) = norm(A,inf).* norm(inv(A),inf);          % A condition number k(A)
    q_norm = norm(q,2);                               % q norm_2 
    v_norm = norm(v,2);                               % v norm_2 
    A_norm_fro = norm(A,'fro');                       % x norm_frobenious

    %---B---%
    y_vec = Ly_b(L,P*v);                              % y from Ly=b
    x_vec = Ux_y(U,y_vec);                            % x from Ux=y
    q_b = x_vec;
    Rel_Error_b(k) = norm(q_b-q',2)./ q_norm;         % Rel error of q'

    %---C---%
    delta_v = v_norm .* 10.^(-3);
    v_c = v + delta_v;
    y_vec = Ly_b(L,P*v_c);
    x_vec = Ux_y(U,y_vec);
    q_c = x_vec;
    Rel_Error_c(k) = norm(q_c-q',2)./ q_norm;

    %---D---%
    delta_A = A_norm_fro .* 10^(-3);
    A_d = delta_A + A;
    [L,U,P] = lu(A_d);
    y_vec = Ly_b(L,P*v);
    x_vec = Ux_y(U,y_vec);
    q_d = x_vec;
    Rel_Error_d(k) = norm(q_d-q',2)./ q_norm;
end

%---Plots---%
figure(1);
subplot(2,1,1);
lg = loglog(h_list,Rel_Error_b,"-",h_list,Rel_Error_c,"-",h_list,Rel_Error_d,"-");
lg(1).LineWidth = 1.5; 
lg(2).LineWidth = 1.5;
lg(3).LineWidth = 1.5;

title('Relative calc error function of h');
ylabel('Relative Error');
xlabel('h(k)');
legend('Relative Error','Relative Error_v','Relative Error_A','Location','northwest');
grid on

subplot(2,1,2);
lg = loglog(h_list,K_A,"-");
lg.LineWidth = 1.5;
title("Condition Number function of h");
ylabel('Cond Num');
xlabel('h(k)');
legend('k(A)[h]');
grid on

movegui(figure(1),"northeast")

%---Functions---%
function A = build_matrix(h,M)
p = 1 ;
A = zeros(M,M);
for m = 1:M
    for n = 1:M
        a_mn = sqrt((h+p*sin(((m*pi)/M))-p*sin(((n*pi)/M))).^2+(p*cos((m*pi)/M)-p*cos((n*pi)/M)).^2);
        A(m,n) = 1./(4*pi*a_mn) ;
    end
end
end

function y = Ly_b(L,b) 
M = length(L);
y = zeros(1,M);
y(1) = b(1) / L(1,1);
sum = 0;
for j=2:M
    for i=1:j-1
        sum = sum + L(j,i) .* y(i);
    end
    y(j) = (b(j) - sum) ./ L(j,j);
    sum = 0;
end
end

function x = Ux_y(U,y)
M = length(U);
x = zeros(1,M);
x(M) = y(M) / U(M,M);
sum = 0;
for i=M-1:-1:1
    for j=M:-1:i+1
        sum = sum + U(i,j).*x(j);
    end
    x(i) = (y(i) - sum) ./ U(i,i);
    sum = 0;
end
end