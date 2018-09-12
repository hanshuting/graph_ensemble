function [L b_ij b_i dpsi_ij dpsi_i] = loss_dispatch(model, psi_ij, psi_i, rho, loss_spec, x)

% provides a convenient way to access all the different loss methods.
% Examples of loss_spec are:
% 'trunc_ul_5'
% 'bbp_ul_1e-5'
% 'pert_ul_1e-5'
% 'trunc_em_5'
% 'em_1e-5'

% 'trunc_trw_ul_5'
% 'trunc_mnf_ul_5'
% 'bbp_trw_ul_1e-5'
% 'bbp_mnf_ul_1e-5'
% 'pert_trw_ul_1e-5'
% 'pert_mnf_ul_1e-5'
% 'trunc_trw/mnf_em_5'
% 'trunc_trw/trw_em_5'
% 'em_trw/mnf_1e-5'
% 'em_trw/trw_1e-5'
% 'pseudo'

% these actually work
% 'trunc_ul_trw_5'
% 

% todo
% 'truncMF' should be an alternative to trunc
% for regular losses it should just substitute mf for trw
% for em it should use meanfield as a lower bound instead of trw
% 'truncTRW' ?!

loss_stuff = regexp(loss_spec,'_','split');

if strcmp(loss_stuff{1},'trunc')
    % need special case for trunc_em
    lossname   = loss_stuff{2};
    inference  = loss_stuff{3};
    maxiter    = str2num(loss_stuff{4});
    convthresh = 0;
    dorec      = 1;
    mode       = 'bprop';
    if strcmp(lossname,'em')
        % call special truncated-em function
        [L b_ij b_i dpsi_ij dpsi_i] = em_trunc(model, psi_ij, psi_i, rho,...
            maxiter, convthresh, x, inference);
    else
        if nargout < 4
            [L b_ij b_i] = implicit_loss(model, psi_ij, psi_i,...
                rho, maxiter, convthresh, x, lossname, dorec, mode, inference);
        else
        [L b_ij b_i dpsi_ij dpsi_i] = implicit_loss(model, psi_ij, psi_i,...
        rho, maxiter, convthresh, x, lossname, dorec, mode, inference);
        end
    end
elseif strcmp(loss_stuff{1},'bbp')
    lossname   = loss_stuff{2};
    inference  = loss_stuff{3};
    maxiter    = 1000; % safety
    convthresh = str2num(loss_stuff{4});
    dorec      = 0;
    mode       = 'bprop';
    if strcmp(lossname,'em')
        error('dont use em with bbp-- just em or trunc_em');
    end
    [L b_ij b_i dpsi_ij dpsi_i] = implicit_loss(model, psi_ij, psi_i,...
    rho, maxiter, convthresh, x, lossname, dorec, mode, inference);
elseif strcmp(loss_stuff{1},'pert') || (length(loss_stuff{1})>4 && strcmp(loss_stuff{1}(1:4),'pert'))
    % can take an argument like
    % 'pert/4' specifying 4-sided differences  or
    % 'pert/4/12.2' specifying 4-sizes differences and emult=12.2
    pertmode  = 2;
    emult     = 1;
    if length(loss_stuff{1})>4
        args = regexp(loss_stuff{1},'/','split');
        pertmode = str2num(args{2});
        if length(args)>2
            emult = str2num(args{3});
        end
    end
    lossname   = loss_stuff{2};
    inference  = loss_stuff{3};
    maxiter    = 1000*10; % safety
    convthresh = str2num(loss_stuff{4});
    dorec      = 0;
    mode       = 'pert';
    if strcmp(lossname,'em')
        error('dont use em with pert-- just em or trunc_em');
    end    
    extra_args = [pertmode emult];

    [L b_ij b_i dpsi_ij dpsi_i] = implicit_loss(model, psi_ij, psi_i,...
    rho, maxiter, convthresh, x, lossname, dorec, mode, inference, extra_args);
elseif strcmp(loss_stuff{1},'em')
    maxiter = 1000; % safety
    inference = loss_stuff{2};
    convthresh = str2num(loss_stuff{3});
    damp = 0;
    [L b_ij b_i dpsi_ij dpsi_i] = em(model, psi_ij, psi_i, rho, maxiter, ...
        convthresh, x, damp, inference);
elseif strcmp(loss_stuff{1},'pseudo')
    %[L b_ij b_i dpsi_ij dpsi_i] = pseudo(model, psi_ij, psi_i, x);
    [L b_ij b_i dpsi_ij dpsi_i] = pseudo_fast(model, psi_ij, psi_i, x);
elseif strcmp(loss_stuff{1},'piecewise')
    [L b_ij b_i dpsi_ij dpsi_i] = piecewise(model, psi_ij, psi_i, x);
elseif strcmp(loss_stuff{1},'piecewise-shared')
    [L b_ij b_i dpsi_ij dpsi_i] = piecewise_shared(model, psi_ij, psi_i, x);
else
    error(['unsupported method: ' loss_stuff{1}])
end