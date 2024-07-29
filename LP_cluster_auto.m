clear all
clc

fileID = fopen('LPresults_cluster2.txt','w');


for n = 4 : 5
    fprintf(fileID,'%2s %1d\n','n =',n);
    fprintf(fileID,'%2s\n','');
    fprintf(fileID,'%1s %1s %6s %6s %6s\n','d       ','k       ','v_small','      v_big','        v_both');
    for d = 2 : 5
        for k = 1 : d-1

            %Operators
            id = eye(d);
            idn = Tensor(id,n);
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

            %Cluster state in the graph basis
            rho = cl*cl';
            rhoG = U*rho*U';

            v1 = LPGMENmixercl(n,d,k,rhoG,0); %Smaller
            v2 = LPGMENmixercl(n,d,k,rhoG,2); %Bigger
            v3 = LPGMENmixercl(n,d,k,rhoG,1); %Both

            v = [d k v1 v2 v3];

            fprintf(fileID,'%1d %8d %13.4f %13.4f %13.4f\n',v);
        end
    end
end

fclose(fileID);