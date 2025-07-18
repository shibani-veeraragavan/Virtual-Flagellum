function Fil = InitialGuess(Fil)
    % Filament.INITIALGUESS()  rearranges linear interpolation as a guess 
    %                          for x^(j+1), i.e.
    %                          x^j = 0.5*( x^(j-1) + x^(j+1) ) .
    %                          Only does position of particle 1, because
    %                          RobotArm does the rest.
        
        % Guess position of particle 1.
        Fil.X(:,1) = 2*Fil.X1(1:3) - Fil.X1(4:6); 
        Utemp = Fil.U;
        Qtemp = Fil.Q;

        for i=1:Fil.N_w
            % Guess Lie algebra elements, and hence quaternions, of all
            % particles.
            Utemp(1:3,i) = 2*Utemp(4:6,i) - Utemp(7:9,i);
            Qtemp(i,1:4) = QuaternionProduct(qexp(Utemp(1:3,i)),...
                                                             Qtemp(i,5:8));
        end
        Fil.U = Utemp;
        Fil.Q = Qtemp;
    end