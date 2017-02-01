function likelihood = compute_sample_likelihood_overcomplete(F,G,logZ,sample,structure)

overcomplete_struct = samples_to_overcomplete(sample, structure);

thetaN = (F*overcomplete_struct.Ut)';
thetaE = (G*overcomplete_struct.Vt)';
linearN = vec(thetaN .* overcomplete_struct.YN);
linearE = vec(thetaE .* overcomplete_struct.YE); 

likelihood = (sum(linearN)) + (sum(linearE)) - logZ;

end