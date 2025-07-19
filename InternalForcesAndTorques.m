function Fil = InternalForcesAndTorques(Fil,dt)
    % Filament.INTERNALFORCESANDTORQUES()  Place gravity, elastic and 
    %                                      constraint forces and torques 
    %                                      on the filament. They are stored
    %                                      in Filament.F

        N = Fil.N_w;  
        Qtemp = Fil.Q;
        Lam = Fil.Lambda;
        dL = Fil.DL;
        MA = Fil.StrainTwist;
        Bend = Fil.KB;
        Twist = Fil.KT;

        % Gravity (in z-direction)
        ForcesAndTorques = [zeros(2,N);...
                            -Fil.weight_per_unit_length*Fil.DL*ones(1,N);...
                            zeros(3,N)]; 
        
        for i=1:Fil.N_w-1
        
            % Constraint forces and torques
            ForcesAndTorques(1:3,i) = ForcesAndTorques(1:3,i) - Lam(:,i);

            ForcesAndTorques(1:3,i+1) = ForcesAndTorques(1:3,i+1) ...
                                      + Lam(:,i);

            t = QuaternionRotation(Qtemp(i,1:4),[1;0;0]);

            ForcesAndTorques(4:6,i) = ForcesAndTorques(4:6,i) - ...
                0.5*dL*[t(2)*Lam(3,i) - t(3)*Lam(2,i);t(3)*Lam(1,i) ...
                                      - t(1)*Lam(3,i);t(1)*Lam(2,i) ...
                                      - t(2)*Lam(1,i)];

            t = QuaternionRotation(Qtemp(i+1,1:4),[1;0;0]);

            ForcesAndTorques(4:6,i+1) = ForcesAndTorques(4:6,i+1) - ...
                0.5*dL*[t(2)*Lam(3,i) - t(3)*Lam(2,i);t(3)*Lam(1,i) ...
                      - t(1)*Lam(3,i);t(1)*Lam(2,i) - t(2)*Lam(1,i)];

            % Elastic and Active Moments
            q = MidpointQ(Qtemp(i,1:4),Qtemp(i+1,1:4));
            dqds = (Qtemp(i+1,1:4) - Qtemp(i,1:4))/dL;
            c = 2*QuaternionProduct([q(1),-q(2),-q(3),-q(4)],dqds);
            C = c(2:4);
            d1 = QuaternionRotation(q,[1;0;0]);
            d2 = QuaternionRotation(q,[0;1;0]);
            d3 = QuaternionRotation(q,[0;0;1]);
            
            M_el = (Twist*C(1) + MA(1,i))*d1 + (Bend*C(2) + MA(2,i))*d2 + ...
                (Bend*C(3) + MA(3,i))*d3;

            ForcesAndTorques(4:6,i) = ForcesAndTorques(4:6,i) + M_el;
            ForcesAndTorques(4:6,i+1) = ForcesAndTorques(4:6,i+1) - M_el;
            
        end
        
        % Update force and moment balances
        Fil.F(:,:) = ForcesAndTorques(:,:);
        
    end
