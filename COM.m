function [x,y,z] = COM(R,M)
%COM function for finding center of mass of multiple body system
%   R is N-by-3 matrix of positions corresponding to N bodies
%   M is N-by-1 column vector of masses corresponding to N bodies (constant
%   mass assumed)
X = 0; Y = 0; Z = 0;
for n = 1:length(M)
    X = X + R(n,1)*M(n);
    Y = Y + R(n,2)*M(n);
    Z = Z + R(n,3)*M(n);
end
x = X/sum(M);
y = Y/sum(M);
z = Z/sum(M);

end

