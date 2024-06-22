%---Q3---%

M = 18;
p = 1 ;
h_val = [1/5,1/2,2,5,10];
h_list = h_val.*((pi * p) ./ M);
q = [3,1,9,0,9,5,8,3,2,3,1,6,3,8,9,5,8,4]'; %charges vector of id's
det_A = zeros(5,1);
Rel_err_q = zeros(5,1);
for t=1:5

    %---A+B---%
     h = h_val(t).*((pi * p) ./ M) ; 
     A = build_matrix(h,M);
     det_A(t) = abs(det(A));
     v = A * q;
     trans_A = transpose(A);
     approx_q = inv(trans_A * A) * trans_A * v;
     Rel_err_q(t) = (norm(approx_q - q)) / norm(q); 
end

%---Plots---%
figure(3);
subplot(2,1,1);
lg = loglog(h_list,det_A,"*-");
lg.LineWidth = 1.5;
title('Least Squares Solution');
xlabel('h');
ylabel('det(A)')
legend('det(A)','Location','southwest');
grid on;

subplot(2,1,2);
lg = loglog(h_list,Rel_err_q,"*-");
lg(1).LineWidth = 1.5;
title('Least Squares Solution');
xlabel('h');
ylabel('rel error')
legend('Rel error of approx q and q','Location','northwest');
grid on;

movegui(figure(3),"southeast")

%---Functions---%
function A = build_matrix(h,M)
p = 1 ;
A = zeros(M,M);
for m = 1:M
    for n = 1:M
        r_mn = sqrt((h+p*sin((m*pi)/M)-p*sin((n*pi)/M)).^2+(p*cos((m*pi)/M)-p*cos((n*pi)/M)).^2);
        A(m,n) = 1 ./ (4*pi*r_mn) ;
    end
end
end