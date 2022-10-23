function C = gram_schmidt(C)
% Orthogonalize C: Gram-Schmidt
use_QR_factorization = 1;
% The QR factorization expresses an m-by-n matrix C as C = Q*R. 
% Q is an m-by-m unitary matrix. 
% R is an m-by-n upper triangular matrix. 
% If the components of C are real numbers, then Q is an orthogonal matrix.
[m,n] = size(C); Q =zeros(m,n); R = zeros(n,n);

if ~use_QR_factorization          % 1. Straight forward method
   for j = 1:n
       v = C(:,j);
       for i = 1:j-1
           R(i,j)=Q(:,i)'*C(:,j);
           v=v-R(i,j)*Q(:,i);     % subtract the projection
       end
       R(j,j) =norm(v);           % v is now perpendicular to all of q1,...q(j-1)
       Q(:,j)=v/R(j,j);
   end
else                              % 2. Using QR Factorization:
   [Q,~] = qr(C);
   Q = Q(:,1:n);
end
C = Q;
end