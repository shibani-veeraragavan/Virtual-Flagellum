function Fil = RobotArm(Fil)
    % Filament.ROBOTARM()  Starting at X(1), construct a beam like a 
    %                      robot arm, by adding a segment a distance 
    %                      Filament.DL  away, with the angle 
    %                      dictated as the average of the two segments' 
    %                      tangents. The process is repeated like this.
    %                      You must call  Filament.InitialSetup  before
    %                      calling this.

        Xtemp = Fil.X;
        Qtemp = Fil.Q;
        dL = Fil.DL;
        for i = 2:Fil.N_w
            Xtemp(:,i) = Xtemp(:,i-1) + ...
                dL/2*(QuaternionRotation(Qtemp(i-1,1:4),[1;0;0]) + ...
                QuaternionRotation(Qtemp(i,1:4),[1;0;0]));
        end
        Fil.X = Xtemp;   
    end