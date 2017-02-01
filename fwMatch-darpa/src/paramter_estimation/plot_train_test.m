function plot_train_test(fig_id,xtrain,ytrain,ytest,type,parameter,measure)

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
        legend(['Log Train ' measure],['Log Test ' measure]);
        ylabel(['Log ' measure]);
    else
        legend(['Train ' measure],['Test ' measure]);
        ylabel(measure);
    end
    grid on;
end