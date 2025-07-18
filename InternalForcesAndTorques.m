function Fil = InternalForcesAndTorques(Fil,dt)
    % Filament.INTERNALFORCESANDTORQUES()  Place gravity, elastic and 
    %                                      constraint forces and torques 
    %                                      on the filament. They are stored
    %                                      in Filament.F

        N = Fil.N_w;  
        Qtemp = Fil.Q;
        Lam = Fil.Lambda;
        dL = Fil.DL;
        ST = Fil.StrainTwist;
        Bend = Fil.KB;
        Twist = Fil.KT;
        NB = Fil.NB;
        NT = Fil.NT;
%         A = Fil.A;
%         B = Fil.B;
%         ST2 = [ST,[0;0;0]];
        % Gravity (in z-direction)
        ForcesAndTorques = [zeros(2,N);...
                            -Fil.weight_per_unit_length*Fil.DL*ones(1,N);...
                            zeros(3,N)]; 
        
        % Elastic torques on segment 1
%         ForcesAndTorques(4:6,1) = ForcesAndTorques(4:6,1) + ...
%                                      Fil.KB/(Fil.N_w * Fil.DL) * [1; 1; 1];

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

            % Elastic, Active and Frictional torques
            
            % Curvature and material frame at current time step, point i+0.5
            q = MidpointQ(Qtemp(i,1:4),Qtemp(i+1,1:4));
            dqds = (Qtemp(i+1,1:4) - Qtemp(i,1:4))/dL;
            c = 2*QuaternionProduct([q(1),-q(2),-q(3),-q(4)],dqds);
            C = c(2:4);
            d1 = QuaternionRotation(q,[1;0;0]);
            d2 = QuaternionRotation(q,[0;1;0]);
            d3 = QuaternionRotation(q,[0;0;1]);
            % Curvature and material frame at previous timestep, point i+0.5
            q = MidpointQ(Qtemp(i,5:8),Qtemp(i+1,5:8));
            dqds = (Qtemp(i+1,5:8) - Qtemp(i,5:8))/dL;
            c = 2*QuaternionProduct([q(1),-q(2),-q(3),-q(4)],dqds);
            C0 = c(2:4);
            d10 = QuaternionRotation(q,[1;0;0]);
            d20 = QuaternionRotation(q,[0;1;0]);
            d30 = QuaternionRotation(q,[0;0;1]);
            % Curvature and material frame two timesteps ago, point i+0.5
            q = MidpointQ(Qtemp(i,9:12),Qtemp(i+1,9:12));
            dqds = (Qtemp(i+1,9:12) - Qtemp(i,9:12))/dL;
            c = 2*QuaternionProduct([q(1),-q(2),-q(3),-q(4)],dqds);
            C00 = c(2:4);
            d100 = QuaternionRotation(q,[1;0;0]);
            d200 = QuaternionRotation(q,[0;1;0]);
            d300 = QuaternionRotation(q,[0;0;1]);
            
            % Elastic and Active Moments
            M_el = (Twist*C(1) + ST(1,i))*d1 + (Bend*C(2) + ST(2,i))*d2 + ...
                (Bend*C(3) + ST(3,i))*d3;
                
%             M_el = QuaternionRotation(q,[Twist*C(1) + ST(1,i); ...
%                 Bend*C(2) + ST(2,i); Bend*C(3) + ST(3,i)]);

            ForcesAndTorques(4:6,i) = ForcesAndTorques(4:6,i) + M_el;
            ForcesAndTorques(4:6,i+1) = ForcesAndTorques(4:6,i+1) - M_el;
            
            % Internal Frictional Moment
            M_di = NT*((3*C(1) - 4*C0(1) + C00(1))*d1/(2*dt) + (3*d1 - 4*d10 + d100)*C(1)/(2*dt)) + ...
                NB*((3*C(2) - 4*C0(2) + C00(2))*d2/(2*dt) + (3*d2 - 4*d20 + d200)*C(2)/(2*dt)) + ...
                NB*((3*C(3) - 4*C0(3) + C00(3))*d3/(2*dt) + (3*d3 - 4*d30 + d300)*C(3)/(2*dt));
            
            ForcesAndTorques(4:6,i) = ForcesAndTorques(4:6,i) + M_di;
            ForcesAndTorques(4:6,i+1) = ForcesAndTorques(4:6,i+1) - M_di;
            
            % Curvature-Controlled Active Moment
%             M_a = dL/2*(A*(C(1)*d1 + C(2)*d2 + C(3)*d3) + B*((3*C(1) - 4*C0(1) + C00(1))*d1/(2*dt) + (3*d1 - 4*d10 + d100)*C(1)/(2*dt) ...
%                 + (3*C(2) - 4*C0(2) + C00(2))*d2/(2*dt) + (3*d2 - 4*d20 + d200)*C(2)/(2*dt) + (3*C(3) - 4*C0(3) + C00(3))*d3/(2*dt) + (3*d3 - 4*d30 + d300)*C(3)/(2*dt)));
%             ST2(:,i) = ST2(:,i) + [sum(M_a.*d1); sum(M_a.*d2); sum(M_a.*d3)];
%             ST2(:,i+1) = ST2(:,i+1) + ST2(:,i) + [sum(M_a.*d1); sum(M_a.*d2); sum(M_a.*d3)];           
            
%             ForcesAndTorques(1:3,i) = ForcesAndTorques(1:3,i) + M_a;
%             ForcesAndTorques(1:3,i+1) = ForcesAndTorques(1:3,i+1) + M_a;
            
        end
        
        % Add curvature-control active moment to StrainTwist vector
%         Fil.StrainTwistJ(:,:) = ST2(:,1:end-1);
        % Update force and moment balances
        Fil.F(:,:) = ForcesAndTorques(:,:);
    end