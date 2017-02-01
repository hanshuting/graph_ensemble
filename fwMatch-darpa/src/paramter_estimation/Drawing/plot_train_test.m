function plot_train_test(fig_id,xtrain,ytrain,ytest,type,parameter,measure,opt)

    figure(fig_id),clf,
    if strcmp(type,'loglog')
        h1 = loglog(xtrain,ytrain);
        hold on,
        h2 = loglog(xtrain,ytest);
    elseif strcmp(type,'semilogx')
        h1 = semilogx(xtrain,ytrain);
        hold on,
        h2 = semilogx(xtrain,ytest);
    elseif strcmp(type,'semilogy')
        h1 = semilogy(xtrain,ytrain);
        hold on,
        h2 = semilogy(xtrain,ytest);
    elseif strcmp(type,'linear')
        h1 = plot(xtrain,ytrain);
        hold on,
        h2 = plot(xtrain,ytest);
    end
    
    set(h1,'LineWidth',2); set(h2,'LineWidth',2);
    set(h1,'Color','b'); set(h2,'Color','r');
    title('Paramter Estimation Performance with \lambda');
    
    if strcmp(type,'semilogx') || strcmp(type,'loglog')
        ylabel(['Log ' parameter]);
    else
        ylabel(parameter);
    end
    
    if strcmp(type,'semilogy') || strcmp(type,'loglog')
        legend(['Log Approx. Test' measure],['Log True Test ' measure]);
        ylabel(['Log ' measure]);
    else
        legend(['Approx. Test ' measure],['True Test ' measure]);
        ylabel(measure);
    end
    grid on;
    
     %% Show optimal lambda
%     if strcmp(opt,'max')
%         [ML,MLI] = max(ytest);
%     elseif strcmp(opt,'min')
%         [ML,MLI] = min(ytest);
%     end
%     opt_xtrain = xtrain(MLI);
%     hold on,
%     if strcmp(type,'semilogx') || strcmp(type,'loglog')
%         xa = [log(opt_xtrain) log(opt_xtrain)];
%     else
%         xa = [opt_xtrain opt_xtrain];
%     end
%     xa = [opt_xtrain opt_xtrain];
%     if strcmp(type,'semilogy') || strcmp(type,'loglog')
%         ya = [log(ML) min(log(yest))];
%     else
%         ya = [ML min(ytest)];
%     end
%     [xaf,yaf] = ds2nfu(xa,ya);
%     xaf, yaf
%     annotation('line',xaf,[yaf(1) 0]);
    
    
end