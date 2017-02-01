function plotHst(h, xFinal, hst)

    figure(h);
    
    iters = [hst.iter];
    
    subplot(2,2,1);
    hist(xFinal);
    title('Final Beliefs');
    
    subplot(2,2,2);
    plot(iters, [hst.fval]);
%     legend('Objective', 'Duality Gap');
    title('Objective and Abs Duality Gap vs. Iter');
    
    subplot(2,2,3);
    semilogy(iters, [hst.trainErr], iters, [hst.relDg]);
    legend('Training Error', 'Rel Duality Gap');
    title('Errors and Rel. Duality Gap vs. Iter');
    
    subplot(2,2,4);
    plot(iters, [hst.stepSz]);
    title('Step size vs. iter');

end

