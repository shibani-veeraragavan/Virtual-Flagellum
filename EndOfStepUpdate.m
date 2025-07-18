function Fil = EndOfStepUpdate(Fil)
    % Filament.ENDOFSTEPUPDATE()  Step in time. Sets historical values of 
    %                             positions (etc.) to their new values.
        
        Fil.X1(4:6) = Fil.X1(1:3);
        Fil.X1(1:3) = Fil.X(:,1);

        Fil.Xm2 = Fil.Xm1;
        Fil.Xm1 = Fil.X;

        Utemp = Fil.U;
        Qtemp = Fil.Q;
        Qtemp(:,9:12) = Qtemp(:,5:8);
        Utemp(7:9,:) = Utemp(4:6,:);
        Utemp(4:6,:) = Utemp(1:3,:);
        
        for i=1:Fil.N_w
            Qtemp(i,5:8) = QuaternionProduct(qexp(Utemp(1:3,i)),...
                                                             Qtemp(i,5:8));
        end

        Fil.U = Utemp;
        Fil.Q = Qtemp;
    end