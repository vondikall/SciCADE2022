function results=logData(A, epsilon, solver)
% LOGDATA Runs all index combinations of A with solver.
% Returns table containing solver info and results.
%   log=logData(A, epsilon, solver)
arguments
    A {mustBeNonempty}
    epsilon {mustBeNonempty}
    solver {mustBeNonempty}
end
n = size(A,1); s = size(A,3);
variable_names_types = [["solver", "string"]; ...
    ["index", "double"]; ...
    ["P","double"]
    ["P_eig_min", "double"]; ...
    ["P_eig_max", "double"]; ...
    ["P_eig_fact", "double"]; ...
    ["Eq_eig_max", "double"]; ...
    ["Solver_info", "string"]; ...
    ["Test_results", "string"]];

results = table('Size',[0,size(variable_names_types,1)],...
    'VariableNames', variable_names_types(:,1),...
    'VariableTypes', variable_names_types(:,2));
for numberOfEquations=1:s
    eqIndex = nchoosek(1:s,numberOfEquations);
    for j=1:size(eqIndex,1)

        % Setup and run solver
        P=sdpvar(n,n);
        constr = [P - epsilon*eye(n) >= 0];
        for i = eqIndex(j,:)
            constr = [constr, A(:,:,i)'*P+P*A(:,:,i)+epsilon*eye(n)<=0];
        end

        opt = sdpsettings('solver', solver, 'warning', 0, 'verbose', 0);
        diagn = optimize(constr,[],opt);

        % Results
        P=double(P); eig_P=eig(P); ii=1;
        for i=eqIndex(j,:)
            eqEigs(ii:(ii+n-1))=eig(A(:,:,i)'*P+P*A(:,:,i));
            ii=ii+n;
        end

        solver_info='';
        if contains(diagn.info,'Successfully')
            solver_info='Success';
            if or(max(eqEigs)>=0,min(eig_P)<=0)
                test_result="Failed test";
            else
                test_result="Passed test";
            end
        elseif contains(diagn.info,'Numerical')
            solver_info='Numerical issues';
        elseif contains(diagn.info,'Infeasible')
            solver_info='Infeasible';
        elseif contains(diagn.info,'Lack of progress')
            solver_info='Lack of progress';
        else
            solver_info='Other';
        end
        if ~contains(diagn.info,'Successfully')
            test_result="Solver failed";
        end

        tmp={solver,{eqIndex(j,:)},P,min(eig_P),max(eig_P),max(eig_P)/min(eig_P),max(eqEigs),solver_info,test_result};
        results=[results;tmp];
    end
end
