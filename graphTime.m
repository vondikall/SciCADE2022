function graphTime(results, saveFigs)
% GRAPHTIME Plots various representations of the data in results. saveFigs
% is optional and automaticly false.
% If saveFig is true it saves four 600p .png figures.
%   graphTime(results)
%   graphTime(results, saveFigs)
arguments
    results {mustBeNonempty}
    saveFigs (1,1) {mustBeNumericOrLogical} = 0
end
solvers = unique(results.solver(:));

figure Name Solver_results
X=categorical(results.Solver_info);
pie(X)

figure Name Test_results
X=categorical(results(results.Solver_info=="Success",:).Test_results);
pie(X)

figure Name Result_overview
result_summary=groupsummary(results,{'solver','Solver_info','Test_results'});
summary=zeros(length(solvers),5);
result_values = ["Failed test", "Passed test", "Numerical issues", "Infeasible", "Lack of progress"];
for i=1:length(solvers)
    right_solver=result_summary.solver==solvers(i);
    for j=1:5
        if j<=2
            tmp=result_summary(and(right_solver, ...
                result_summary.Test_results==result_values(j)),:).GroupCount;
        else
            tmp=result_summary(and(right_solver, ...
                result_summary.Solver_info==result_values(j)),:).GroupCount;
        end
        if isempty(tmp)
            summary(i,j)=0;
        else
            summary(i,j)=tmp;
        end

    end
end
X = categorical(solvers);
X = reordercats(X,solvers);
b = bar(X,summary);
for i=1:5   
    xtips = b(i).XEndPoints;
    ytips = b(i).YEndPoints;
    labels = string(b(i).YData);
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
end
legend(result_values,'Location','eastoutside')


figure Name Eigenvalue_info
results_passed = results(results.Test_results=="Passed test",:);
for i=1:length(solvers)
    P_fact=results_passed.P_eig_fact(results_passed.solver==solvers(i));
    Eq_eig=results_passed.Eq_eig_max(results_passed.solver==solvers(i));
    scatter(P_fact,Eq_eig,'filled'); hold on;
end
xlabel('\lambda_{P,max}/\lambda_{P,min}')
ylabel('\lambda_{Eq,max}')
set(gca,'xscale','log','yscale','log')
legend(solvers, 'Location','northwest')

if saveFigs
    figure(1)
    print('solver_results_all','-dpng', '-r600')
    figure(2)
    print('test_results_all','-dpng', '-r600')
    figure(3)
    print('Results_bar','-dpng', '-r600')
    figure(4)
    print('Eigenvalue_info_all','-dpng', '-r600')
end