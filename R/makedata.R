# Generation of simulations

# contains anomalies and mean changes, synced across dimensions
simData = function(P, N, s_s, t0_s, t1_s, sps = FALSE, times = NULL){
  # Y will be an P x N matrix
  # s_s is the number of true signals with changes
  # t0_s is the number of time points with anomalies
  # t1_s is the number of time points with mean changes
  # if sps, not all dimensions change together
  # times are the change locations
  
  K_s = min(P, N)
  
  # dimensions with change points
  sigs = sample(1:K_s, s_s)
  
  # change points are at least 1 apart
  if (is.null(times)){
    times = sample(seq(2, N-1, by = 2), t0_s + t1_s)    
  }else{
    times = sample(times, length(times))
  }
  times0 = times1 = numeric(0)
  if (t0_s > 0){
    times0 = sort(times[1:t0_s])
  }
  times1 = sort(times[! times %in% times0])
  times = sort(times)
  
  M_s = matrix(runif(P*K_s, min = -1, max = 1), nrow = P, ncol = K_s)
  M_s[,-sigs] = 0
  S_s = matrix(numeric(N*K_s), nrow = N, ncol = K_s)
  for (s in K_s){
    S_s[,s] = runif(1, 1, 5)*sample(c(-1, 1), 1)
  }

  if (length(times) > 0){
    for (i in 1:length(times)){
      ch = (times[i] %in% times1)
      start_t = times[i]
      end_t = ifelse(i == length(times), N+1, times[i+1])
      if (ch==1){
        # mean change
        if (sps){
          # sample dimensions that change, no change if FALSE
          ch_TF = numeric(s_s)
          while (sum(ch_TF) == 0){
            ch_TF = sample(c(TRUE, FALSE), s_s, replace = TRUE)  
          }
          ch_sigs = sigs[ch_TF]
        }else{
          ch_sigs = sigs
        }
        for (s in ch_sigs){
          S_s[start_t:N, s] = S_s[start_t-1,s] + 
            rep(runif(1, 1, 5)*sample(c(-1, 1), 1), N+1-start_t)
        }
      }else if (ch==0){
        # anomaly
        if (sps){
          # sample dimensions that change, no change if FALSE
          ch_TF = numeric(s_s)
          while (sum(ch_TF) == 0){
            ch_TF = sample(c(TRUE, FALSE), s_s, replace = TRUE)  
          }
          ch_sigs = sigs[ch_TF]
        }else{
          ch_sigs = sigs
        }
        for (s in ch_sigs){
          S_s[start_t:(end_t-1),s] = c(S_s[start_t-1,s] + runif(1, 1, 5)*sample(c(-1, 1), 1), 
                                       rep(S_s[start_t-1,s], end_t-start_t-1))
        }
      }
    }
  }

  rowvars_s = runif(P, 0.1, 5)
  E_s = matrix(rnorm(P*N, 0, rep(sqrt(rowvars_s), each=N)), 
               nrow = P, ncol = N, byrow=TRUE)
  
  Y = M_s%*%t(S_s) + E_s
  
  return(list(Y = Y, M = M_s, S = S_s, times0 = times0, times1 = times1,   
              T0 = M_s%*%t(S_s), psi = rowvars_s))
}
