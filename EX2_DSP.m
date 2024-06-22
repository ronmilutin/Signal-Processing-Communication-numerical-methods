disp('316389584')
clear all
%% PART A Q4 %%
%4b
fprintf('\nw0=0\n')
ERG1_X(0)
fprintf('\nw0=0.01\n')
ERG1_X(0.01)
fprintf('\nw0=0.1\n')
ERG1_X(0.1)
fprintf('\nw0=pi/6\n')
ERG1_X(pi/6)

%4c
fprintf('\nw0=0\n')
ERG2_X(0)
fprintf('\nw0=pi/2\n')
ERG2_X(pi/2)
fprintf('\nw0=pi\n')
ERG2_X(pi)
%functions%
%a
function  X = GEN_X(w0, N)
A = randn;
B = randn;
n = 1:N;
z_n=randn(1,N);
X=A*cos(w0*n)+B*sin(w0*n)+z_n;
end
%b
function ERG1_X(w0)
rafs = 100;
N=1000;
expectation = 0;
differences= zeros(rafs,1);
for i=1:rafs
    X = GEN_X(w0, N);
    est_expext=mean(X);
    differences(i)=est_expext-expectation;
    end
fprintf ('mean error:%8.5f\n', mean(differences));
fprintf ('Error std: %8.5f\n', std(differences));
end

%c
function ERG2_X(w0)
rafs = 100;
N=1000;
expectation = 0;
k=0:10;
true_rx=cos(w0*k);
true_rx0=2;
differences2=zeros(rafs,11);
for i=1:rafs
    X=GEN_X(w0,N);
    for k=0:10
        est_rx(k+1)=sum(X(1:end-10).*X(1+k:end-10+k))/(N-10);
    end
    differences2(i,:)= est_rx - true_rx;
end

fprintf('k:        ');
fprintf('%8d',0:10);
fprintf('\n');
fprintf('mean error= ');
fprintf('%8.5f', mean(differences2));
fprintf('\n');
fprintf('Error std:');
fprintf('%8.5f', std(differences2));
fprintf('/n');
end



