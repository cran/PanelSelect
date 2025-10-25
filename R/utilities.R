# a few utility functions
panel.select.env = new.env(parent = emptyenv())
panel.select.env$iter = 1

resetIter = function(){
    panel.select.env$iter = 1
}

addIter = function(){
    panel.select.env$iter = panel.select.env$iter + 1
}

.onUnload <- function (libpath) {
    library.dynam.unload("PanelSelect", libpath)
}

getVarSE = function(res, gvar=NULL, verbose=0){
    # Greene: For small and moderate sized samples, Hessian is preferable
    # Simulation shows that the original asymmetric hessian produces much smaller SE and the 95% coverage probability is closer to 95%

    if(!is.null(gvar)){
        res$var_bhhh = gvar$var
        rownames(res$var_bhhh) = rownames(res$hessian)
        colnames(res$var_bhhh) = colnames(res$hessian)
        res$se_bhhh = sqrt(diag(res$var_bhhh))
    }

    res$var = tryCatch( solve(-res$hessian), error = function(e){
        cat('Hessian not invertible: ', e$message, '\n')
        m = NA * diag(length(res$par))
        rownames(m) = colnames(m) = names(res$par)
        m
    } )

    res$variance_type = 'Hessian'
    if(any(is.na(res$var)) || any(diag(res$var)<0)) {
        if(!is.null(res$var_bhhh)){
            if(all(!is.na(res$var_bhhh)) && all(diag(res$var_bhhh)>0)) {
                cat('Negative/Missing values in Hessian-style variance matrix, using BHHH instead\n')
                if(verbose>-1) print(res$var)
                res$variance_type = 'BHHH'
                res$var = res$var_bhhh
            } else {
                warning('Negative/Missing values in both Hessian-style and BHHH-style variance matrices')
                if(verbose>-1) {
                    cat('Hessian and BHHH variance matrices below:\n')
                    print(res$var)
                    print(res$var_bhhh)
                }
            }
        }else warning('Negative values in Hessian-style variance matrix, no BHHH available\n')
    }
    res$se = sqrt(diag(res$var))
    res$gtHg = matrix(res$gradient, nrow=1) %*% res$var %*% matrix(res$gradient, ncol=1)
    res$AIC = -2*res$LL + 2 * length(res$par)
    res$BIC = -2*res$LL + log(res$n_obs) * length(res$par)
    res$n_par = length(res$par)

    res
}

compile = function(par, se){
    lci = par - qnorm(0.975)*se
    uci = par + qnorm(0.975)*se
    z = par/se
    p = 1 - pchisq(z^2, 1)
    cbind(estimate=par,se=se,z=z,p=p,lci=lci,uci=uci)
}

transCompile = function(res, trans_vars=NULL, trans_types=NULL){
    if(length(trans_vars)!=length(trans_types))
        stop('Numbers of variables and functions are different for transformation')

    # Retain the raw estimate and se
    par = res$par
    se = res$se
    var = res$var
    npar = length(par)

    lci = par - qnorm(0.975)*se
    uci = par + qnorm(0.975)*se

    if(length(trans_vars)>=1){
        for(i in 1:length(trans_vars)){
            var_name = trans_vars[i]
            type = trans_types[i]
            g = 1 # gradient of transformation function

            # Note: delta method is not reliable when the original parameter is not normal
            # https://www.statisticshowto.com/delta-method-definition/#:~:text=Disadvantages,LePage%20%26%20Billard%2C%201992).
            if(type=='exp'){ # >0
                par[var_name] = exp(par[var_name])
                lci[var_name] = exp(lci[var_name])
                uci[var_name] = exp(uci[var_name])
                g = par[var_name]
            }else if(type=='correlation'){ # [-1, 1]
                tau = par[var_name]
                par[var_name] = 1 - 2/(exp(tau)+1)
                lci[var_name] = 1 - 2/(exp(lci[var_name])+1)
                uci[var_name] = 1 - 2/(exp(uci[var_name])+1)
                g = 2*exp(tau)/(exp(tau)+1)^2
            }else if(type=='logit'){ # [0,1]
                q_eta = par[var_name]
                par[var_name] = plogis(q_eta)
                lci[var_name] = plogis(lci[var_name])
                uci[var_name] = plogis(uci[var_name])
                g = dlogis(q_eta)
            }else{
                cat(sprintf('No transformation is done for parameter %s with type %s\n', var_name, type))
            }
            # update SE and variance
            se[var_name] = g * se[var_name]
            # var(theta) = t(dtheta/dpar) %*% var(par) %*% dtheta/dpar
            # dtheta/dpar is a diagonal matrix with g on the diagonal
            var[var_name, ] = g * var[var_name, ]
            var[, var_name] = g * var[, var_name]

            ix = which(names(par)==var_name)
            names(se)[ix] = names(par)[ix] = names(lci)[ix] = names(uci)[ix] = names(trans_vars)[i]
            colnames(var)[ix] = rownames(var)[ix] = names(se)[ix]
        }
    }
    z = par/se
    p = 1 - pchisq(z^2, 1)
    res$estimates = cbind(estimate=par,se=se,z=z,p=p,lci=lci,uci=uci)
    res$par = par
    res$se = se
    res$var = var
    res
}

dbinorm = function(x1, x2, rho){
    z = -(x1^2 + x2^2 - 2*rho*x1*x2) / (2*(1-rho^2))
    exp(z) / (2*pi*sqrt(1-rho^2))
}

