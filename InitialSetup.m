function Fil = InitialSetup(Fil,FirstSegPos)
    % Filament.INITIALSETUP(FirstSegPos,StrainTwist,B,Sp,omega_b,...
    %                                            weight_per_unit_length,mu)
    %     Set first segment at position FirstSegPos, then use RobotArm to 
    %     construct the rest of the segments. Properties of the filament 
    %     are as the names of the inputs suggest.
        
        Fil.X(:,1) = FirstSegPos;
        Fil.X1 = [FirstSegPos;FirstSegPos];
        N = Fil.N_w;
        % The orientation quaternions are initialised to the identity by 
        % default.
        Fil.Q(1:N,1) = ones(N,1);
        Fil.Q(1:N,5) = ones(N,1);
        Fil.Q(1:N,9) = ones(N,1);
        Fil = RobotArm(Fil);
        Fil.Xm1 = Fil.X;
        Fil.Xm2 = Fil.X;
    end
