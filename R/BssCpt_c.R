# R wrapper

bsscpt = function(Y, K = 5, burn_in = 500, numiter = 3000, a.psi = 1, b.psi = 1){
  # stage 1 with partial model
  output1 = 'error'
  attempt = 1
  while (class(output1) == 'character' && attempt <= 5){
    output1 = tryCatch({
      CptC(Y = Y, K = K, d = 1, numiter = numiter, a_psi = a.psi, b_psi = b.psi)          
    }, error = function(e){
      'error'
    })
    attempt = attempt + 1
  }
  # stage 2 with full model
  if (class(output1) != 'character'){
    output1ests = CptCEsts(output1, burn_in = burn_in)
    init = setInit_sps(output1ests)    
    output2 = 'error'
    attempt = 1
    while (class(output2) == 'character' && attempt <= 5){
      output2 = tryCatch({
        CptCinit(Y = Y, K = K, numiter = numiter, a_psi = a.psi, b_psi = b.psi, 
                 init = init, useInit = TRUE) 
      }, error = function(e){
        'error'
      })
      attempt = attempt + 1
    }
    if (class(output2) != 'character'){
      output2ests = CptCinitEsts(output2, burn_in = burn_in)  
      res = list(Mest = output2ests[['Mest']], V0est = output2ests[['V0est']], V1est = output2ests[['V1est']],
                 Sest = output2ests[['Sest']], psiest = output2ests[['psiest']],
                 g0est = output2ests[['g0est']], g1est = output2ests[['g1est']],
                 cpt0 = output2ests[['cpt0']], cpt1 = output2ests[['cpt1']])
    }else{
      res = list(Mest = init[['M']], V0est = init[['V0']], V1est = init[['V1']],
                 Sest = init[['S']], psiest = init[['psi']],
                 cpt0 = init[['cpt0']], cpt1 = init[['cpt1']])
    }
  }else{
    print('Error in gibbs sampling')
    res = list(cpt0 = integer(0), cpt1 = integer(0))
  }

  return(res)
}

# helper functions

CptCEsts = function(output, burn_in){
  V = output[['V']]
  M = output[['M']]
  S = output[['S']]
  tau = output[['tau']]
  xi = output[['xi']]
  lam = output[['lam']]
  v = output[['v']]
  phi = output[['phi']]
  w = output[['w']]
  gamma = output[['gamma']]
  alpha = output[['alpha']]
  psi = output[['psi']]
  
  numiter = dim(M)[3]
  samples = (burn_in + 1):numiter
  N = dim(V)[1]
  theta = array(dim=c(dim(M)[1],dim(S)[1],length(samples)))
  g = array(dim=c(dim(phi)[1],length(samples)))
  for(i in 1:length(samples)){
    theta[,,i] = M[,,samples[i]]%*%t(S[,,samples[i]])
    idx = apply(abs(V[,,i]),1,which.max)
    g_i = c()
    for(j in 1:N){
      g_i = c(g_i,V[j,idx[j],i])
    }
    g[,i] = g_i
  }
  thetaest = apply(theta,1:2,median)
  gest = apply(g,1,median)
  
  Mest = apply(M[,,burn_in:numiter],1:2,median)
  Sest = apply(S[,,burn_in:numiter],1:2,median)
  Vest = apply(V[,,burn_in:numiter],1:2,median)
  
  tauest = median(tau[burn_in:numiter])
  xiest = median(xi[burn_in:numiter])
  lamest = apply(lam[,burn_in:numiter],1,median)
  vest = apply(v[,burn_in:numiter],1,median)
  phiest = apply(phi[,burn_in:numiter],1,median)
  west = apply(w[,burn_in:numiter],1,median)
  gammaest = apply(gamma[,,burn_in:numiter],1:2,median)
  alphaest = apply(alpha[,,burn_in:numiter],1:2,median)
  psiest = apply(psi[,burn_in:numiter],1,median)
  
  return(list(tauest = tauest, xiest = xiest, lamest = lamest, vest= vest, 
              phiest = phiest, west = west, gammaest = gammaest, alphaest = alphaest,
              Mest = Mest, Vest = Vest, Sest = Sest, psiest = psiest, 
              thetaest = thetaest, gest = gest))
  
}

CptCinitEsts = function(output, burn_in){
  V0 = output[['V0']]
  V1 = output[['V1']]
  M = output[['M']]
  S = output[['S']]
  tau0 = output[['tau0']]
  tau1 = output[['tau1']]
  lam0 = output[['lam0']]
  lam1 = output[['lam1']]
  phi0 = output[['phi0']]
  phi1 = output[['phi1']]
  gamma0 = output[['gamma0']]
  gamma1 = output[['gamma1']]
  psi = output[['psi']]
  
  numiter = dim(M)[3]
  samples = (burn_in + 1):numiter
  N = dim(V0)[1]
  theta = array(dim=c(dim(M)[1],dim(S)[1],length(samples)))
  g0 = array(dim=c(dim(phi0)[1],length(samples)))
  g1 = array(dim=c(dim(phi1)[1],length(samples)))
  for(i in 1:length(samples)){
    theta[,,i] = M[,,samples[i]]%*%t(S[,,samples[i]])
    idx0 = apply(abs(V0[,,i]),1,which.max)
    idx1 = apply(abs(V1[,,i]),1,which.max)
    g0_i = c()
    g1_i = c()
    for(j in 1:N){
      g0_i = c(g0_i,V0[j,idx0[j],i])
      g1_i = c(g1_i,V1[j,idx1[j],i])
    }
    g0[,i] = g0_i
    g1[,i] = g1_i
  }
  thetaest = apply(theta,1:2,median)
  g0est = abs(apply(g0, 1, median))
  g1est = abs(apply(g1, 1, median))
  
  Mest = apply(M[,,burn_in:numiter],1:2,median)
  Sest = apply(S[,,burn_in:numiter],1:2,median)
  V0est = apply(V0[,,burn_in:numiter],1:2,median)
  V1est = apply(V1[,,burn_in:numiter],1:2,median)
  
  tau0est = median(tau0[burn_in:numiter])
  tau1est = median(tau1[burn_in:numiter])
  lam0est = apply(lam0[,burn_in:numiter],1,median)
  lam1est = apply(lam1[,burn_in:numiter],1,median)
  phi0est = apply(phi0[,burn_in:numiter],1,median)
  phi1est = apply(phi1[,burn_in:numiter],1,median)
  gamma0est = apply(gamma0[,,burn_in:numiter],1:2,median)
  gamma1est = apply(gamma1[,,burn_in:numiter],1:2,median)
  psiest = apply(psi[,burn_in:numiter],1,median)
  
  #### Select change points by kde
  
  if (max(g0est) <= 1){
    cutoff0 = max(g0est)
  }else{
    cutoff0 = max(g0est)
    den = density(g0est, kernel = 'rectangular')
    den_diff = diff(den[['y']])
    pos = 1
    for (j in 1:length(den_diff)){
      if (pos == 1 && den_diff[j] < 0 && den[['x']][j] > 0){
        pos = 2
      }
      if (pos == 2 && den_diff[j] >= 0 && den[['y']][j] < 1e-10 &&
          den[['x']][j] > 0){
        cutoff0 = den[['x']][j]
        break
      }
    }  
  }
  if (max(g1est) <= 1){
    cutoff1 = max(g1est)
  }else{
    cutoff1 = max(g1est)
    den = density(g1est, kernel = 'rectangular')
    den_diff = diff(den[['y']])
    pos = 1
    for (j in 1:length(den_diff)){
      if (pos == 1 && den_diff[j] < 0 && den[['x']][j] > 0){
        pos = 2
      }
      if (pos == 2 && den_diff[j] >= 0 && den[['y']][j] < 1e-10 &&
          den[['x']][j] > 0){
        cutoff1 = den[['x']][j]
        break
      }
    }  
  }
  cpt0 = (1:N)[g0est > cutoff0]
  cpt1 = (1:N)[g1est > cutoff1]
  # remove 1 since that's the starting value
  cpt0 = cpt0[! cpt0 %in% 1]
  cpt1 = cpt1[! cpt1 %in% 1]
  
  return(list(tau0est = tau0est, tau1est = tau1est, lam0est = lam0est, lam1est = lam1est, 
              phi0est = phi0est, phi1est = phi1est, gamma0est = gamma0est, gamm1est = gamma1est, 
              Mest = Mest, V0est = V0est, V1est = V1est, Sest = Sest, psiest = psiest, 
              thetaest = thetaest, g0est = g0est, g1est = g1est,
              cpt0 = cpt0, cpt1 = cpt1))
  
}

# process results from partial model for initialization
# cutoff for zero values of gest determined by kde
setInit_sps = function(res1){
  # res1 are results from partial model
  
  gest = res1[['gest']]
  phiest = res1[['phiest']]
  west = res1[['west']]
  gammaest = res1[['gammaest']]
  alphaest = res1[['alphaest']]
  Vest = res1[['Vest']]
  
  # cutoff for zero values of gest
  gest_abs = abs(gest)
  if (max(gest_abs) <= 1){
    cutoff = max(gest_abs)
  }else{
    cutoff = max(gest_abs)
    den = density(gest_abs, kernel = 'rectangular')
    den_diff = diff(den[['y']])
    pos = 1
    for (j in 1:length(den_diff)){
      if (pos == 1 && den_diff[j] < 0 && den[['x']][j] > 0){
        pos = 2
      }
      if (pos == 2 && den_diff[j] >= 0 && den[['y']][j] < 1e-10 &&
          den[['x']][j] > 0){
        cutoff = den[['x']][j]
        break
      }
    }       
  }
  
  # locating possible change point locations
  index = which(gest_abs > cutoff)
  index_0 = which(gest_abs <= cutoff)
  # index of possible additive outliers
  index_diff = index[which(diff(index) == 1)]
  index_anom = integer(0)
  i = 1
  while (i <= length(index_diff)){
    pn = sign(gest[index_diff[i]])*sign(gest[index_diff[i]+1])
    if (pn == -1){
      anom = index_diff[i]
      inseq = TRUE      
    }
    i = i + 1
    if (pn == -1){
      while (inseq){
        if (index_diff[i] - index_diff[i-1] == 1 && 
            sign(gest[index_diff[i]])*sign(gest[index_diff[i]+1]) == -1 &&
            i < length(index_diff)){
          anom = c(anom, index_diff[i])
          i = i + 1
        }else{
          inseq = FALSE
        }
      }
      index_anom = c(index_anom, seq(anom[1], anom[length(anom)], by = 2))      
    }
  }
  index_conseq = c()
  for (x in index_anom){
    index_conseq = c(index_conseq, x, x + 1)  
  }
  if (!is.null(index_conseq)){
    index_conseq = sort(unique(index_conseq))
  }
  # index of possible level shifts
  index_mc = index[! index %in% index_conseq]
  
  # remove 1 since that's the starting value
  cpt0 = index_anom[! index_anom %in% 1]
  cpt1 = index_mc[! index_mc %in% 1]
  
  # phi
  N = length(phiest)
  med = median(phiest[index_0])
  phi0est = phi1est = rep(med, N)
  phi0est[index_anom] = phiest[index_anom]
  phi1est[index_mc] = phiest[index_mc]
  
  # w
  med = median(west[index_0])
  w0est = w1est = rep(med, N)
  w0est[index_anom] = west[index_anom]
  w1est[index_mc] = west[index_mc]
  
  # gamma
  K = ncol(gammaest)
  med = median(gammaest[index_0,])
  gamma0est = gamma1est = matrix(rep(med, N*K), nrow = N, ncol = K)
  gamma0est[index_anom,] = gammaest[index_anom,]
  gamma1est[index_mc,] = gammaest[index_mc,]
  
  # alpha
  med = median(alphaest[index_0,])
  alpha0est = alpha1est = matrix(rep(med, N*K), nrow = N, ncol = K)
  alpha0est[index_anom,] = alphaest[index_anom,]
  alpha1est[index_mc,] = alphaest[index_mc,]
  
  # V
  V0est = V1est = matrix(numeric(N*K), nrow = N, ncol = K)
  V0est[index_anom,] = Vest[index_anom,]
  V1est[index_mc,] = Vest[index_mc,]
  
  
  return(list(tau = res1[['tauest']], xi = res1[['xiest']], lam = res1[['lamest']], v = res1[['vest']], 
              phi0 = phi0est, phi1 = phi1est, w0 = w0est, w1 = w1est,
              gamma0 = gamma0est, gamma1 = gamma1est, alpha0 = alpha0est, alpha1 = alpha1est,
              M = res1[['Mest']], V0 = V0est, V1 = V1est, S = res1[['Sest']], psi = res1[['psiest']],
              cpt0 = cpt0, cpt1 = cpt1))  
}


# dynamic programming to prune estimated change points
# by sum of squared errors from mean
prettycpt = function(S, cpt0, cpt1){

  # find points from cpt1
  N = ncol(S)
  Q = length(cpt1)
  if (Q <= 2){
    cpt1_new = cpt1
  }else{
    endloc = c(cpt1[2:Q]-1, N)
    sse_total = matrix(NA, Q, Q)
    sse_right = matrix(NA, Q, Q) # saves sse on right of change point
    cpt1_prev = matrix(NA, Q, Q) # saves previous change location
    # iterate through each number of changes
    for (q in 1:Q){
      if (q == 1){
        # iterate through each ending location
        for (n in 1:Q){
          # iterate through each candidate change location
          best_sse = Inf
          for (l in 0:(n-1)){ # l here is col position in sse_total
            cpt1_l = cpt1[l+1]
            leftseg = S[, 1:(cpt1_l-1), drop = FALSE]
            rightseg = S[, cpt1_l:endloc[n], drop = FALSE]
            leftsse = sum((leftseg - rowMeans(leftseg))**2)
            rightsse = sum((rightseg - rowMeans(rightseg))**2)
            # fill sse_right
            if (l != 0){
              sse_right[l, n] = rightsse
            }
            if ((leftsse + rightsse) < best_sse){
              best_sse = leftsse + rightsse
              best_l = cpt1_l
            }
          }
          sse_total[q, n] = best_sse
          cpt1_prev[q, n] = best_l
        }      
      }else{
        # q > 1 
        # iterate through each ending location 
        for (n in q:Q){
          # iterate through each candidate change location
          best_sse = Inf
          for (l in (q-1):(n-1)){ 
            totalsse = sse_total[q-1, l] + sse_right[l, n]
            if (totalsse < best_sse){
              best_sse = totalsse
              best_l = l
            }          
          }
          sse_total[q, n] = best_sse
          cpt1_prev[q, n] = best_l
        }
      }
    }
    # pick change points with mean difference of sse as threshold
    thresh = mean(diff(sse_total[,Q]))
    qest = max(which(diff(sse_total[,Q]) <= thresh)) + 1
    colnum = cpt1_prev[qest, Q]
    cpt1_new = endloc[colnum] + 1
    if (qest > 2){
      for (q in (qest-1):2){
        colnum = cpt1_prev[q, colnum]
        cpt1_new = c(endloc[colnum]+1, cpt1_new)
      }      
    }
    cpt1_new = c(cpt1_prev[1, colnum], cpt1_new)
  }
  
  # find points from cpt0
  sse0 = c()
  endloc_new = c(0, cpt1_new-1, N)
  for (i in 1:(length(endloc_new)-1)){
    startseg = endloc_new[i] + 1
    endseg = endloc_new[i+1]
    cpt0_seg = cpt0[cpt0 >= startseg & cpt0 <= endseg]
    if (length(cpt0_seg) > 0){
      segmean = rowMeans(S[, startseg:endseg, drop = FALSE])
      for (c0 in cpt0_seg){
        sse0 = c(sse0, sum((S[, c0] - segmean)**2))
      }
    }
  }
  # pick change points with mean difference of sse0 as threshold
  cpt0_new = cpt0[sse0 >= mean(sse0)]

  return(list(cpt0 = cpt0_new, cpt1 = cpt1_new))
}
