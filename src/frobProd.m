function y = frobInnerProduct(A, B)
% frobInnerProduct  Frobenius inner product of A and B
%   http://stackoverflow.com/questions/8031628/octave-matlab-efficient-calc-of-frobenius-inner-product

	y = A(:).'*B(:);

end

