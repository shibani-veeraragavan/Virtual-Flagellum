function out = epsilon(i,j,k)
out = (i*j*k==6)*(j-i)*(k-i)*(k-j)/2;
end