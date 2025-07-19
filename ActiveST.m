function Fil = ActiveST(Fil, omega_a, k_a, A, nt)
% Calculate sinusoidal active moment wave vector, and place into Filament.StrainTwist 
        
        L = Fil.Length;
        KB = Fil.KB;
        N = Fil.N_w;
        dL = Fil.DL;
        ST = zeros(3,N-1);
        s = dL:dL:L; % Arclength
        for i=1:N-1
            ST(1,i) = A(1)*KB/L*cos(k_a(1)*s(i) - omega_a(1)*nt);
            ST(2,i) = A(2)*KB/L*cos(k_a(2)*s(i) - omega_a(2)*nt);
            ST(3,i) = A(3)*KB/L*cos(k_a(3)*s(i) - omega_a(3)*nt);
        end
        Fil.StrainTwist = ST;
    end
