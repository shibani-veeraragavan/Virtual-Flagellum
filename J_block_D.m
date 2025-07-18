function [D] = J_block_D(Filament,dt,mu)
% J_BLOCK_D(Filament,mu)   Provides the terms relating the velocity
%                          constraints to the Lagrange multipliers.
% Rows correspond to m=1...N-1, Columns correspond to n=1...N-1.
% Input: Filament (object), dt (scalar), mu (scalar).
% Output: C_mn = df^3_m/dF_n-0.5 (Size 3(N-1) x 3(N-1)) 

N_w = Filament.N_w;

D = zeros(3*(N_w-1));
fac = zeros(N_w);

for m=N_w:-1:1
    
    fac(m) = 1/(6*pi*mu*Filament.R(m));
    % fac(m) = -dt/(9*pi*mu*Filament.R(m));
    
end

% At every row block I, add -fac_m*[I] at n=m+1; [I] is the identity tensor.
% m = I+1;
for I=1:N_w-2
    
    i = 3*(I-1); % Ending point of (m-1) block of rows.
    j = i+3;
    
    D(i+1,j+1) = -fac(I+1);
    D(i+2,j+2) = -fac(I+1);
    D(i+3,j+3) = -fac(I+1);
    
    % Adds +fac_m*[I] to all rows at n=1 (columns 1-3) for row blocks
    % m=1...N-2.
    % Delete next 3 lines %
    D(i+1,1) = fac(1);
    D(i+2,2) = fac(1);
    D(i+3,3) = fac(1);
    
end
% Add +fac_n*[I] to row block m=N-1 at n=2 (columns 1-3).
% Delete next 3 lines %
D(end-2,1) = fac(1);
D(end-1,2) = fac(1);
D(end,3) = fac(1);

% Add +fac_n*[I] at all diagonal elements (n=m).
for i=1:3*(N_w-1)
    
    I = 1 + ceil(i/3); 
    
    D(i,i) = D(i,i) + fac(I);
    
end

end