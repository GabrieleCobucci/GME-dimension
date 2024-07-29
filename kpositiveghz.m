function L = kpositiveghz(F,k,n,d,cut,perm)

%-------------------------------------------------------------------------%
%This function applies the k-positive generalised reduction map:
%(Id_S \otimes \Lambda_{\bar{S}})[F_{S|\bar{S}}],
%with \Lambda[X] = \tr[X]Id - 1/k*X

%Inputs:
% - F: state density matrix (diagonal in the GHZ basis)
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

%GHZ state
ghz = 0;
for a = 0 : d-1
    ghz = ghz + 1/sqrt(d)*Tensor(comp{a+1},n);
end

%Unitary from GHZ basis (computational) to graph basis (diagonal)
U = 0;
for j = 0 : d^n-1
    L = toSeveralBases(j,d*ones(1,n));
    oper = Z^(L(1));
    for t = 2 : length(L)
        oper = Tensor(oper,X^(L(t)));
    end
    ghzj = oper*ghz;
    U = U + idn(:,j+1)*ghzj';
end

%Generalised reduction map
L = U*PermuteSystems(Tensor(idA,PartialTrace(U'*F*U,cut,d*ones(1,n))),perm)*U' - 1/k*F;

end