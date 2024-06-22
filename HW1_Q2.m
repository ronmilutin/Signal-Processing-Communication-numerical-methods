%---Q2---%

M = 18;
p = 1;
h_val = [5,2,1];
q = [3,1,9,0,9,5,8,3,2,3,1,6,3,8,9,5,8,4]';
for t = 1:3

    %---A+B---%
    h = (pi * p)./ (h_val(t).*M);
    A = build_matrix(h,M);
    v = A * q;
    if t == 1
        [Rel_dist_A,Rel_err_A,iter_A] = Gauss_Seidel(A,v,q); %h=pi/5M
        %---C---%
        [Rel_dist_Jac,Rel_err_Jac,iter_Jac] = Jacobi(A,v,q); %h=pi/5M
        %---D---%
        A_d = build_matrix_d(h,M);
        v_d = A_d * q;
        [Rel_dist_D,Rel_err_D,iter_D] = Jacobi(A_d,v_d,q); %h=pi/5M
    end
    if t == 2
        [Rel_dist_B,Rel_err_B,iter_B] = Gauss_Seidel(A,v,q); %h=pi/2M
    end
    if t == 3 
        [Rel_dist_C,Rel_err_C,iter_C] = Gauss_Seidel(A,v,q); %h=pi/M
    end          
end

%---Plots---%
figure(2);
subplot(4,2,[1,2])
plot_A = semilogy(iter_A,Rel_dist_A,iter_A,Rel_err_A,'-*');
plot_A(1).LineWidth = 1.5;
plot_A(2).LineWidth = 1.5;
msg2 = sprintf('Took %d iterations',length(iter_A));
msg2_1 = sprintf('Final Rel error of q^k and q is %d',Rel_err_A(end));
text(2,10^-4,msg2,'Color','green','FontSize',12);
text(2,10^-4.5,msg2_1,'Color','green','FontSize',12);
title(' Gauss-Seidel For h=pi/5M');
xlabel('Iterations');
legend('Rel Destination of q^k and q^(k-1)','Rel Real Error of q^k and q^(k-1)');

subplot(4,2,[3,4])
plot_B = semilogy(iter_B,Rel_dist_B,iter_B,Rel_err_B,'-*');
plot_B(1).LineWidth = 1.5;
plot_B(2).LineWidth = 1.5;
msg3 = sprintf('Took %d iterations',length(iter_B));
msg3_1 = sprintf('Final Rel error of q^k and q is %d',Rel_err_B(end));
text(2,10^-4,msg3,'Color','green','FontSize',12);
text(2,10^-4.5,msg3_1,'Color','green','FontSize',12);
title('Gauss-Seidel For h=pi/2M');
legend('Rel Destination of q^k and q^(k-1)','Rel Real Error of q^k and q^(k-1)');

subplot(4,2,[5,6])
plot_C = semilogy(iter_C,Rel_dist_C,iter_C,Rel_err_C,'-*');
plot_C(1).LineWidth = 1.5;
plot_C(2).LineWidth = 1.5;
msg4 = sprintf('Solution Does not Converge');
text(2,10^-4,msg4,'Color','green','FontSize',12);
title('Gauss-Seidel For h=pi/M');
legend('Rel Destination of q^k and q^(k-1)','Rel Real Error of q^k and q^(k-1)');

subplot(4,2,7)
plot_Jac = semilogy(iter_Jac,Rel_dist_Jac,iter_Jac,Rel_err_Jac,'-*');
plot_Jac(1).LineWidth = 1.5;
plot_Jac(2).LineWidth = 1.5;
msg5 = sprintf('Solution Does not Converge');
text(10,10^6,msg5,'Color','green','FontSize',12);
title('Jacobi for h=pi/5M');
legend('Rel Destination of q^k and q^(k-1)','Rel Real Error of q^k and q^(k-1)');

subplot(4,2,8)
plot_D = semilogy(iter_D,Rel_dist_D,iter_D,Rel_err_D,'-*');
plot_D(1).LineWidth = 1.5;
plot_D(2).LineWidth = 1.5;
msg6 = sprintf('Took %d Iterations',length(iter_D));
msg6_1 = sprintf('Final Relative error of q^k and q is %d',Rel_err_D(end));
text(1.5,10^-4.5,msg6,'Color','green','FontSize',12);
text(1.5,10^-4.7,msg6_1,'Color','green','FontSize',12);
title('Jacobi for h=pi/5M with new matrix A');
legend('Rel Destination of q^k and q^(k-1)','Rel Real Error of q^k and q^(k-1)');

movegui(figure(2),"northwest")

%---Functions---%
function A = build_matrix(h,M)
p = 1 ;
A = zeros(M,M);
for m = 1:M
    for n = 1:M
        a_mn = sqrt((h+p*sin((m*pi)/M)-p*sin((n*pi)/M)).^2+(p*cos((m*pi)/M)-p*cos((n*pi)/M)).^2);
        A(m,n) = 1 ./ (4*pi*a_mn);
    end
end
end

function A = build_matrix_d(h,M)
p = 1 ;
A = zeros(M,M);
for m = 1:M
    for n = 1:M
        a_mn = (h+p*sin((m*pi)/M)-p*sin((n*pi)/M)).^2+(p*cos((m*pi)/M)-p*cos((n*pi)/M)).^2;
        A(m,n) = 1 ./ (4*pi*a_mn) ;
    end
end
end

function [Rel_dist,Rel_err,iter_plot,q_k] = Gauss_Seidel(A,v,q)
    L = tril(A,-1);
    D = diag(diag(A)); 
    Q = L + D ;
    neg_u = Q - A; 
    inv_Q = inv(Q);
    G = inv_Q * neg_u ;         % -u * (L+D)^(-1)
    C =  inv_Q * v;             % (L+D)^(-1) * v
    Err_Endurance = 10 ^ (-3);
    q_k = C;                    % for q^(1), k=1
    Rel_dist = zeros;
    Rel_err = zeros;
    iter = 1;   
    iter_plot = zeros;
    err = max(abs(q-q_k));
    max_iter = 500;             %Limit Iterations
    while abs(err) > Err_Endurance && iter <=max_iter  
       q_k_minus_1 = q_k;
       q_k = G*(q_k_minus_1) + C; % q^(k)=-u*(L+D)^(-1) * q^(k-1) + (L+D)^(-1)*v
       err = norm(q-q_k,'inf');
       Rel_dist(iter) = norm(q_k - q_k_minus_1, 'inf') / norm(q_k_minus_1, 'inf');
       Rel_err(iter) = norm(q_k - q,'inf') / norm(q,'inf');
       iter_plot(iter) = iter;
       iter = iter + 1;
    end
end

function [Rel_dist,Rel_err_real,iter_plot,q_k] = Jacobi(A,v,q)   
    D = diag(diag(A));
    Q = D ;
    I = eye(18);
    inv_Q = inv(Q);
    G = I - (inv_Q * A);
    C = inv_Q * v;
    Err_Endurance = 10 ^ (-3);
    q_k = C;                    %for q^(1), k=1
    Rel_dist = zeros;
    Rel_err_real = zeros;
    iter = 1;   
    iter_plot = zeros;
    err = max(abs(q-q_k));
    max_iter = 400;             %Limit Iterations
    while abs(err) > Err_Endurance && iter <=max_iter  
       q_k_minus_1 = q_k;
       q_k = G*(q_k_minus_1) + C; % q^(k)=-u*(L+D)^(-1) * q^(k-1) + (L+D)^(-1)*v
       err = norm(q-q_k,'inf');
       Rel_dist(iter) = norm(q_k - q_k_minus_1, 'inf') / norm(q_k_minus_1, 'inf');
       Rel_err_real(iter) = norm(q_k - q,'inf') / norm(q,'inf');
       iter_plot(iter) = iter;
       iter = iter + 1;
    end
end