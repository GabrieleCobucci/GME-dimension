function v = LPGMENmixercl(n,d,k,phi,option)

%-------------------------------------------------------------------------%
%This function evaluates the linear program (LP) for a target pure
%n-partite high-dimensional linear cluster state by using YALMIP.

%Inputs:
% - n: number of parties
% - d: dimension of local subsystems
% - k: GME-dimension
% - phi: pure target density matrix
% - option: application mode of the generalised reduction map
% -- 0: Partial trace on the smallest subset
% -- 2: Partial trace on the biggest subset
% -- 1: Map applied on both subsets

%Output:
% - v: maximum visibility at which a noisy GHZ state has GME-dimension = k
%-------------------------------------------------------------------------%

%Identity
id = eye(d);

%Number of bipartitions
s = Stirling2nd(n,2);

%Set of bipartitions
bipartitions = SetPartition(n,2);

%Constraints
C = [];

%Visibility
v = sdpvar(1);

%Total state and diagonal variables
sigma = 0;
for m = 1 : s
    tau{m} = diag(sdpvar(d^n,1));
    sigma = sigma + tau{m};
    C = [C, diag(tau{m}) >= 0];
end

%Reduction map constraint
for j = 1 : s
    %Discriminate between smallest and biggest subsets
    b = cell2mat(bipartitions{j}(1));
    c = cell2mat(bipartitions{j}(2));
    if length(b) <= length(c)
        set = [b];
        comp = [c];
    else
        set = [c];
        comp = [b];
    end
    %Map on the smallest subset
    if option == 0
        order1 = [set comp];
        [~,perm1] = sort(order1);
        rm = kpositivecluster(tau{j},k,n,d,set,perm1);
        C = [C, rm >= 0];
    %Map on the both subsets
    elseif option == 1
        order1 = [set comp];
        [~,perm1] = sort(order1);
        order2 = [comp set];
        [~,perm2] = sort(order2);
        rm1 = kpositivecluster(tau{j},k,n,d,set,perm1);
        rm2 = kpositivecluster(tau{j},k,n,d,comp,perm2);
        C = [C, rm1 >= 0, rm2 >= 0];
    %Map on the biggest subset
    elseif option == 2
        order1 = [comp set];
        [~,perm1] = sort(order1);
        rm = kpositivecluster(tau{j},k,n,d,comp,perm1);
        C = [C, rm >= 0];
    else
        error('Optimization option not valid: choose a value in {0,1,2}');
    end
end

%Target state constraint
C = [C, sigma == v*phi + (1-v)*Tensor(id,n)/d^n];


%SolveSDP
disp('Options')
ops=sdpsettings('solver','scs', 'cachesolvers', 1);
diagnostic=solvesdp(C,-v,ops)

v=double(v)

end