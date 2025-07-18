function [F] = J_block_F(Filament,dt,mu)
% J_BLOCK_F(Filament,dt,mu)   Provides the terms relating the Lie algebra 
%                             update equations to the Lagrange multipliers.
% Rows correspond to m=1...N, Columns correspond to n=2...N.
% Input: Filament (object), dt (scalar), mu (scalar).
% Output: C_mn = df^2_m/dF_n-0.5 (Size 3N x 3(N-1))

N_w = Filament.N_w;

DL = Filament.DL;

R = Filament.R;

U = Filament.U;

Q = Filament.Q;

F = zeros(3*N_w,3*(N_w-1));

for J=1:N_w-1 % Column number, J = n-1; Equivalent to n=2:N_w
    % At every n, we will have entries corresponding to m=n-1 and m=n.
    pos = 3*(J-1);
    
    % For m=n-1
    fac = dt*DL/(24*pi*mu*R(J)^3);
    % fac = -dt*DL/(24*pi*mu*R(J)^3); %
    Dmat = rcross(QuaternionRotation(Q(J,1:4),[1;0;0]))'; % [D1_m x]
    umat = rcross(U(1:3,J))'; % [u_m x]
    
    F(pos+1:pos+3,pos+1:pos+3) = fac*(Dmat - 0.5*umat*Dmat + umat^2*Dmat/12); 
    
    % For m=n
    fac = dt*DL/(24*pi*mu*R(J+1)^3);
    % fac = -dt*DL/(24*pi*mu*R(J+1)^3); %
    Dmat = rcross(QuaternionRotation(Q(J+1,1:4),[1;0;0]))'; 
    umat = rcross(U(1:3,J+1))';
    
    F(pos+4:pos+6,pos+1:pos+3) = fac*(Dmat - 0.5*umat*Dmat + umat^2*Dmat/12); 
    
end

end
