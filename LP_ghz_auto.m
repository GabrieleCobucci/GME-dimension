clear all
clc

fileID = fopen('LPresults_ghz.txt','w');

for n = 5 : 7
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

            %GHZ state in the graph basis
            rho = ghz*ghz';
            rhoG = U*rho*U';

            v1 = LPGMENmixerghz(n,d,k,rhoG,0); %Smaller
            v2 = LPGMENmixerghz(n,d,k,rhoG,2); %Bigger
            v3 = LPGMENmixerghz(n,d,k,rhoG,1); %Both

            v = [d k v1 v2 v3];

            fprintf(fileID,'%1d %8d %13.4f %13.4f %13.4f\n',v);
        end
    end
end

fclose(fileID);