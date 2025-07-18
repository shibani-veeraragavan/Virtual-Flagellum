function [C] = J_block_C(Filament,fac)
% J_BLOCK_C(Filament,fac)  Provides the terms related to the derivative of
%                          the velocity constraints with respect to the Lie
%                          algebra elements.                        
% Rows correspond to m=2...N, Columns correspond to n=1...N.
% Input: Filament (object), fac (scalar, = dL/2).
% Output: C_mn = df^3_m/du_n (Size 3(N-1) x 3N) 

N_w = Filament.N_w;

Q = Filament.Q;

C = zeros(3*(N_w-1),3*N_w);

% First row, I=1 (m=2) has 3x3 blocks of [xD1_n] at n=1,2.
% D1_n = QuaternionRotation(Q_n^t+1, e_x). [xD1_n] = rcross(D1_n).

C(1:3,1:6) = fac*[rcross(QuaternionRotation(Q(1,1:4),[1;0;0])),...
                  rcross(QuaternionRotation(Q(2,1:4),[1;0;0]))];

for I=2:N_w-1   %Rows for m=3...N
    
    m = I+1; 
    i = 3*(I-1)+ 1;     % Starting point of row I;
    
    % Bring forward the contents of previous row, which occupy C(I-1,
    % 1:I) to current row C(I,1:I).
    C(i:i+2,1:i+2) = C(i-3:i-1,1:i+2);
    % Add 2 new blocks at C(I,I) and C(I, I+1) corresponding to [xD1_m-1]
    % and [xD1_m].
    C(i:i+2,i:i+5) = C(i:i+2,i:i+5) ...
                 + fac*[rcross(QuaternionRotation(Q(m-1,1:4),[1;0;0])), ...
                        rcross(QuaternionRotation(Q(m,1:4),[1;0;0]))];
end

end
