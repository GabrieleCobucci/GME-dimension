function L = kpositivecluster(F,k,n,d,cut,perm)

%-------------------------------------------------------------------------%
%This function applies the k-positive generalised reduction map:
%(Id_S \otimes \Lambda_{\bar{S}})[F_{S|\bar{S}}],
%with \Lambda[X] = \tr[X]Id - 1/k*X

%Inputs:
% - F: state density matrix (diagonal in the cluster state basis)
% - k: GME-dimension
% - n: number of parties
% - d: dimension of local subsystems
% - cut: set on which apply the k-positive operator \Lambda
% - perm: permutation vector to restore the order of the parties after the
% partial trace

%Output:
% - L: operator after the application of the map
%-------------------------------------------------------------------------%

%Identity
id = eye(d);

%Identity of the subset S
a = length(cut);
idA = eye(d^a);

%Global identity
idn = Tensor(id,n);

%Generalized Pauli operators
X = GenPauli(1,0,d);
Z = GenPauli(0,1,d);

%Computational basis
for l = 0 : d-1
    comp{l+1} = id(:,l+1);
end

%n-qudit cluster state
cl = 0;
for l = 0 : d^n-1
    L = toSeveralBases(l,d*ones(1,n));
    term = comp{L(1)+1};
    for m = 2 : n
        term = Tensor(term,Z^(L(m-1))*comp{L(m)+1});
    end
    cl = cl + (1/d)^(n/2)*term;
end

%Unitary from n-qudit cluster (computational) basis to graph diagonal basis
U = 0;
for j = 0 : d^n-1
    L = toSeveralBases(j,d*ones(1,n));
    oper = 1;
    for m = 1 : n-2
        oper = Tensor(oper,Z^(L(m)));
    end
    for m = n-1 : n
        oper = Tensor(oper,X^(L(m)));
    end
    graphj = oper*cl;
    U = U + idn(:,j+1)*graphj';
end

%Generalised reduction map
L = U*PermuteSystems(Tensor(idA,PartialTrace(U'*F*U,cut,d*ones(1,n))),perm)*U' - 1/k*F;

end