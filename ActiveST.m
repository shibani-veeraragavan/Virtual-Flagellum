function Fil = ActiveST(Fil, omega_a, k_a, phi, A, nt)
    % Filament.ActiveMoment(omega_b, k_b, m_0, nt)  Calculate sinusoidal 
    %               active moment wave vector, and place into Filament.mA 
        
        L = Fil.Length;
        KB = Fil.KB;
        N = Fil.N_w;
        dL = Fil.DL;
        ST = zeros(3,N-1);
        s = dL:dL:L; % Arclength
        for i=1:N-1
            if A(4) == 0
                ST(1,i) = A(1)*KB/L*cos(k_a(1)*s(i) - omega_a(1)*nt) + KB/L*phi(1);
                ST(2,i) = A(2)*KB/L*cos(k_a(2)*s(i) - omega_a(2)*nt) + KB/L*phi(2);
                ST(3,i) = A(3)*KB/L*cos(k_a(3)*s(i) - omega_a(3)*nt) + KB/L*phi(3);
            elseif A(4) == 1
                ST(1,i) = A(1)*KB/L*s(i)*cos(k_a(1)*s(i) - omega_a(1)*nt) + KB/L*phi(1);
                ST(2,i) = A(2)*KB/L*s(i)*cos(k_a(2)*s(i) - omega_a(2)*nt) + KB/L*phi(2);
                ST(3,i) = A(3)*KB/L*s(i)*cos(k_a(3)*s(i) - omega_a(3)*nt) + KB/L*phi(3);
            elseif A(4) == 2
                ST(3,i) = A(3)*KB/L*k_a(3)^2*sin(k_a(3)*s(i) - 2*pi*nt)/sqrt(1 - (A(3)*KB/L)^2*k_a(3)^2*cos(k_a(3)*s(i) - 2*pi*nt)^2) + phi(3);
            end
        end
        Fil.StrainTwist = ST;
    end