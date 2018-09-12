function  plotpath(filename)

coeffs  = mmread(filename);
lambdas = mmread([filename,'_lambda']);

coeffs  = coeffs(any(coeffs,2),:);

semilogx(lambdas, coeffs');

