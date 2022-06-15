function main(options)
% ------------------------- SciCADE2022 LMI code -------------------------
% Written by Stefania Andersen Aradottir, saa20 (at) hi.is as a
% supplementary to the poster 'Shortfalls of semidefinite programming when
% applied to real life mechanical systems'
%     main(options)
%     options:
%       saveFig         [1|0 (0)] Saves all figures as 600p .png
%       load_results    [1|0 (1)] If false it will run YALMIP with
%                       predefined solvers
%       close_and_clear [1|0 (0)] If true will run close all and clc before
%                       running code
%
% Additional notes: When not loading results you need to have YALMIP set up
% with sedumi, mosek and sdpt3
% Requires that the following files be saved in the same folder or
% accessible via a saved path:
%   ACC90_benchmark_system.mat
%   results.mat                     (if load_results = true)
%   graphTime.m
%   logData.m                       (if load_results = false)
arguments
    options.saveFig (1,1) {mustBeNumericOrLogical} = 0
    options.load_results (1,1) {mustBeNumericOrLogical} = 1
    options.close_and_clear (1,1) {mustBeNumericOrLogical} = 0
end

if options.close_and_clear
    close all; clc;
end

% Config
saveFig=options.saveFig;
load_results=options.load_results;


% Loading data
load('ACC90_benchmark_system.mat')
[s, n] = size(L); epsilon=1e-10;
A = zeros(n,n,s);
A_eig=zeros(n,s);
A_legend={};

for i=1:s
    A(:,:,i) = [A_org-L(i,1:4)'*C_org, N_Matrix; -L(i,5)*C_org,0];
    A_eig(:,i) = eig(A(:,:,i));
    A_legend{i}=sprintf("A_{%i}",i);
end

% Ploting system eigenvalues.
figure Name System_eigenvalues
subplot(2,1,1); hold on;
subplot(2,1,2); hold on;
for i=1:s
    if i<=s/2
        subplot(2,1,1)
        scatter(real(A_eig(:,i)),imag(A_eig(:,i)), 'filled')
    else
        subplot(2,1,2)
        scatter(real(A_eig(:,i)),imag(A_eig(:,i)), 'filled')
    end
end
subplot(2,1,1)
xline(0,':','imaginary','LabelHorizontalAlignment','left')
yline(0,':','real','LabelHorizontalAlignment','left')
legend(A_legend(1:s/2),'Location','bestoutside');
subplot(2,1,2)
xline(0,':','imaginary','LabelHorizontalAlignment','left')
yline(0,':','real','LabelHorizontalAlignment','left')
legend(A_legend((s/2+1):s),'Location','bestoutside');
if saveFig
    print('A_eig','-dpng', '-r600')
end

if ~load_results
    fprintf(['You chose to run the solvers. ', ...
        'This might take up to 10 minuest\n'])
    fprintf(['NB: Code assumes you have YALMIP set up along ', ...
        'with the following solvers: sedumi, mosek and sdpt3.\n'])
    m=input('Do you want to continue? Y/N [Y]:','s');
    if m=='N'
        load_results=1;
    else
        solvers = ["sedumi", "mosek", "sdpt3"];
        fprintf('Solver started ')
        results = logData(A,epsilon,solvers(1));
        fprintf('-')
        for i=2:length(solvers)
            tmp = logData(A,epsilon,solvers(i));
            results=[results;tmp];
            fprintf('-')
        end
        for i=1:s
            P=lyap(A(:,:,i)',eye(n));
            eig_P=eig(P);
            eqEigs(:)=eig(A(:,:,i)'*P+P*A(:,:,i));
            if or(max(eqEigs)>=0,min(eig_P)<=0)
                test_result="Failed test";
            else
                test_result="Passed test";
            end
            tmp={'lyap',{i},P,min(eig_P),max(eig_P),max(eig_P)/min(eig_P),max(eqEigs),"Success",test_result};
            results=[results;tmp];
        end
    end
    fprintf('\n')
    graphTime(results, saveFig)
else
    load('results.mat')
    graphTime(results, saveFig)
end
m=input('Do you want save the workspace? Y/N [N]:','s');
if m=='Y'
    save('workspace.mat');
end