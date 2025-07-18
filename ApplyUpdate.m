function Fil = ApplyUpdate(Fil,u)
    % Filament.APPLYUPDATE(u)  Alg 2, Line 5. Does  X = X + u  for state
    %                          vector X.
        N = Fil.N_w;
        Fil.X(:,1) = Fil.X(:,1) + u(1:3);
        Utemp = Fil.U;
        Qtemp = Fil.Q;
        Lam = Fil.Lambda;
        for i=1:N
            Utemp(1:3,i) = Utemp(1:3,i) + u(3*i+1:3*(i+1));
            Qtemp(i,1:4) = QuaternionProduct(qexp(Utemp(1:3,i)),...
                                                             Qtemp(i,5:8));
            if i<N
                Lam(:,i) = Lam(:,i) + u(3*(N+i)+1:3*(N+i+1));
            end
        end
        Fil.U = Utemp;
        Fil.Q = Qtemp;
        Fil.Lambda = Lam;
    end