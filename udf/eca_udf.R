eca_two<-function (seriesA, seriesB, alpha = 0.05, delT = 0, sym = FALSE, 
          tau = 0, sigtest = "poisson", reps = 1000) 
{
  if (length(table(seriesA)) > 2) {
    print("|-------------    ERROR #1    -------------|")
    print("|       Time seriesA is not binary         |")
    print("|      (or the number of events=0)!        |")
    print("| Use CC.binarize() to preprocess seriesA  |")
    print("|------------------------------------------|")
    return()
  }
  if (length(table(seriesB)) > 2) {
    print("|-------------    ERROR #1    -------------|")
    print("|       Time seriesB is not binary         |")
    print("|      (or the number of events=0)!        |")
    print("| Use CC.binarize() to preprocess seriesA  |")
    print("|------------------------------------------|")
    return()
  }
  if (tau < 0 || delT < 0) {
    print("|-------------    ERROR #2    -------------|")
    print("|    The offset (tau) or delta T (delT)    |")
    print("|              is negative.                |")
    print("|------------------------------------------|")
    return()
  }
  if (length(seriesA) != length(seriesB)) {
    print("|-------------    ERROR #4    -------------|")
    print("|    lengths of series A and B differ      |")
    print("|------------------------------------------|")
    return()
  }
  if (!is.vector(seriesA) && !is.matrix(seriesA)) {
    print("|-------------    ERROR #5    -------------|")
    print("| series A is neither vector nor matrix    |")
    print("|------------------------------------------|")
    return()
  }
  if (!is.vector(seriesB) && !is.matrix(seriesB)) {
    print("|-------------    ERROR #5    -------------|")
    print("| series B is neither vector nor matrix    |")
    print("|------------------------------------------|")
    return()
  }
  if (any(is.na(seriesA)) && sigtest == "surrogate") {
    print("|-------------    ERROR #6    -------------|")
    print("|          seriesA contains NAs.           |")
    print("| surrogate test not allowed, use poisson  |")
    print("|------------------------------------------|")
    return()
  }
  if (any(is.na(seriesB)) && sigtest == "surrogate") {
    print("|-------------    ERROR #6    -------------|")
    print("|          seriesB contains NAs.           |")
    print("| surrogate test not allowed, use poisson  |")
    print("|------------------------------------------|")
    return()
  }
  if (any(is.na(seriesB)) && sigtest == "shuffle") {
    print("|-------------    ERROR #6    -------------|")
    print("|    seriesB or seriesA contains NAs.      |")
    print("| shuffle test not allowed, use poisson    |")
    print("|------------------------------------------|")
    return()
  }
  if (sum(is.na(seriesA)) == length(seriesA) || sum(is.na(seriesB)) == 
      length(seriesB)) {
    CA_out = list(NA, NA, NA, NA,NA, NA, NA, NA, NA,NA)
    names(CA_out) = c("NH precursor", "NH trigger", "p-value precursor", 
                      "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                      "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")
    
    return(CA_out)
  }
  seriesA[is.na(seriesB)] = NA
  seriesB[is.na(seriesA)] = NA
  bindata = matrix(NA, 2, length(seriesA))
  bindata[1, ] = seriesA
  bindata[2, ] = seriesB
  rownames(bindata) = c("seriesA", "seriesB")
  Tlen = length(bindata[1, !is.na(bindata[1, ])])
  steplen = length(seriesA)
  N_A = as.numeric(Tlen - table(bindata[1, ] == 1)[1])
  N_B = as.numeric(Tlen - table(bindata[2, ] == 1)[1])
  K_prec = 0
  for (step in 1:steplen) {
    if (is.na(bindata[1, step])) {
      next
    }
    if (bindata[1, step] == 1) {
      start = step - tau - (delT)
      if (sym == TRUE) {
        end = step - tau + (delT)
      }
      if (sym == FALSE) {
        end = step - tau
      }
      if (start < 1 & end >= 1) {
        start = 1
      }
      if (start < 1 & end < 1) {
        next
      }
      if (start > steplen) {
        next
      }
      if (end > steplen) {
        end = steplen
      }
      if (is.element(1, bindata[2, start:end])) {
        K_prec = K_prec + 1
      }
    }
  }
  CRprec = K_prec/N_A
  K_trigg = 0
  for (step in 1:steplen) {
    if (is.na(bindata[2, step])) {
      next
    }
    if (bindata[2, step] == 1) {
      if (sym == TRUE) {
        start = step + tau - (delT)
      }
      if (sym == FALSE) {
        start = step + tau
      }
      end = step + tau + (delT)
      if (start < 1 & end >= 1) {
        start = 1
      }
      if (start < 1 & end < 1) {
        next
      }
      if (start > steplen) {
        next
      }
      if (end > steplen) {
        end = steplen
      }
      if (is.element(1, bindata[1, start:end])) {
        K_trigg = K_trigg + 1
      }
    }
  }
  CRtrigg = K_trigg/N_B
  if (sigtest == "poisson") {
    if (sym == FALSE) {
      Pprec = 0
      for (Ktmp in K_prec:N_A) {
        Ptmp = choose(N_A, Ktmp) * ((1 - ((1 - ((delT + 
                                                   1)/(Tlen - tau)))^N_B))^Ktmp) * (((1 - ((delT + 
                                                                                              1)/(Tlen - tau)))^N_B)^(N_A - Ktmp))
        Pprec = Pprec + Ptmp
        Eprec = 1 - (1 - delT/(Tlen - tau))^N_B
      }
      Ptrigg = 0
      for (Ktmp in K_trigg:N_B) {
        Ptmp = choose(N_B, Ktmp) * ((1 - ((1 - ((delT + 
                                                   1)/(Tlen - tau)))^N_A))^Ktmp) * (((1 - ((delT + 
                                                                                              1)/(Tlen - tau)))^N_A)^(N_B - Ktmp))
        Ptrigg = Ptrigg + Ptmp
        Etrigg = 1 - (1 - delT/(Tlen - tau))^N_A
        
      }
      
    }
    if (sym == TRUE) {
      Pprec = 0
      delTsym = delT + delT + 1
      for (Ktmp in K_prec:N_A) {
        Ptmp = choose(N_A, Ktmp) * ((1 - ((1 - (delTsym/(Tlen)))^N_B))^Ktmp) * 
          (((1 - (delTsym/(Tlen)))^N_B)^(N_A - Ktmp))
        Pprec = Pprec + Ptmp
        Eprec = 1 - (1 - delT/(Tlen - tau))^N_B
        
      }
      Ptrigg = 0
      for (Ktmp in K_trigg:N_B) {
        Ptmp = choose(N_B, Ktmp) * ((1 - ((1 - (delTsym/(Tlen)))^N_A))^Ktmp) * 
          (((1 - (delTsym/(Tlen)))^N_A)^(N_B - Ktmp))
        Ptrigg = Ptrigg + Ptmp
        Etrigg = 1 - (1 - delT/(Tlen - tau))^N_A
        
      }
    }
  }
  if (sigtest == "wt.surrogate") {
    if (delT == 0) {
      delT = 1
    }
    seriesA = CC.ts2es(seriesA)
    seriesB = CC.ts2es(seriesB)
    span = (max(seriesA$span[1], seriesB$span[1]):min(seriesA$span[2], 
                                                      seriesB$span[2]))
    Tlen = length(span)
    seriesA_sort = sort(seriesA$es)
    wtA = rep(0, N_A)
    wtA[1] = seriesA_sort[1] - seriesA$span[1]
    for (i in 2:N_A) {
      wtA[i] = seriesA_sort[i] - seriesA_sort[i - 1]
    }
    seriesB_sort = sort(seriesB$es)
    wtB = rep(0, N_B)
    wtB[1] = seriesB_sort[1] - seriesB$span[1]
    for (i in 2:N_B) {
      wtB[i] = seriesB_sort[i] - seriesB_sort[i - 1]
    }
    surdist = matrix(NA, 2, reps)
    for (surno in 1:reps) {
      surA = sample(wtA, size = 1)
      for (i in 2:length(seriesA$es)) {
        tmp = surA[i - 1] + sample(wtA, size = 1)
        if (tmp > span[Tlen]) {
          break
        }
        else {
          surA[i] = tmp
        }
      }
      surB = sample(wtB, size = 1)
      for (i in 2:length(seriesB$es)) {
        tmp = surB[i - 1] + sample(wtB, size = 1)
        if (tmp > span[Tlen]) {
          break
        }
        else {
          surB[i] = tmp
        }
      }
      K_prec_sur = 0
      for (step_a in 1:N_A) {
        if (!is.element(surA[step_a], span)) {
          next
        }
        for (step_b in 1:N_B) {
          if (!is.element(surB[step_b], span)) {
            next
          }
          if (sym == FALSE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(0:delT))) {
              K_prec_sur = K_prec_sur + 1
              break
            }
          }
          if (sym == TRUE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(-delT:delT))) {
              K_prec_sur = K_prec_sur + 1
              break
            }
          }
        }
      }
      surdist[1, surno] = K_prec_sur/N_A
      K_trigg_sur = 0
      for (step_b in 1:N_B) {
        if (!is.element(surB[step_b], span)) {
          next
        }
        for (step_a in 1:N_A) {
          if (!is.element(surA[step_a], span)) {
            next
          }
          if (sym == FALSE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(0:delT))) {
              K_trigg_sur = K_trigg_sur + 1
              break
            }
          }
          if (sym == TRUE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(-delT:delT))) {
              K_trigg_sur = K_trigg_sur + 1
              break
            }
          }
        }
      }
      surdist[2, surno] = K_trigg_sur/N_B
    }
    Pprec = 1 - ecdf(surdist[1, ])(CRprec)
    Ptrigg = 1 - ecdf(surdist[2, ])(CRtrigg)
    Eprec = mean(surdist[1, ])
    Etrigg = mean(surdist[2, ])
  }
  if (sigtest == "shuffle.surrogate") {
    surdist = matrix(NA, 2, reps)
    span = seq(1:steplen)
    for (surno in 1:reps) {
      surA = sample(span, size = N_A)
      surB = sample(span, size = N_B)
      K_prec_sur = 0
      for (step_a in 1:N_A) {
        for (step_b in 1:N_B) {
          if (sym == FALSE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(0:delT))) {
              K_prec_sur = K_prec_sur + 1
              break
            }
          }
          if (sym == TRUE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(-delT:delT))) {
              K_prec_sur = K_prec_sur + 1
              break
            }
          }
        }
      }
      surdist[1, surno] = K_prec_sur/N_A
      K_trigg_sur = 0
      for (step_b in 1:N_B) {
        for (step_a in 1:N_A) {
          if (sym == FALSE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(0:delT))) {
              K_trigg_sur = K_trigg_sur + 1
              break
            }
          }
          if (sym == TRUE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(-delT:delT))) {
              K_trigg_sur = K_trigg_sur + 1
              break
            }
          }
        }
      }
      surdist[2, surno] = K_trigg_sur/N_B
    }
    Pprec = 1 - ecdf(surdist[1, ])(CRprec)
    Ptrigg = 1 - ecdf(surdist[2, ])(CRtrigg)
    Eprec = mean(surdist[1, ])
    Etrigg = mean(surdist[2, ])
    
  }
  sig_testprec = logical
  if (Pprec >= alpha) {
    sig_testprec = TRUE
  }
  else {
    sig_testprec = FALSE
  }
  sig_testtrigg = logical
  if (Ptrigg >= alpha) {
    sig_testtrigg = TRUE
  }
  else {
    sig_testtrigg = FALSE
  }
  CA_out = list(sig_testprec, sig_testtrigg, Pprec, Ptrigg, 
                CRprec, CRtrigg, N_A, N_B, Eprec,Etrigg)
  names(CA_out) = c("NH precursor", "NH trigger", "p-value precursor", 
                    "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                    "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")
  return(CA_out)
}

eca_three <- function (seriesA, seriesB, seriesC, alpha = 0.05, delT = 0, 
                       delT.cond = 0, sym = FALSE, sym.cond = FALSE, tau = 0, tau.cond = 0, 
                       sigtest = "poisson", reps = 1000) 
{
  if (length(table(seriesA)) > 2) {
    print("|-------------    ERROR #1    -------------|")
    print("|       Time seriesA is not binary         |")
    print("|      (or the number of events=0)!        |")
    print("| Use CC.binarize() to preprocess seriesA  |")
    print("|------------------------------------------|")
    return()
  }
  if (length(table(seriesB)) > 2) {
    print("|-------------    ERROR #1    -------------|")
    print("|       Time seriesB is not binary         |")
    print("|      (or the number of events=0)!        |")
    print("| Use CC.binarize() to preprocess seriesB  |")
    print("|------------------------------------------|")
    return()
  }
  if (length(table(seriesC)) > 2) {
    print("|-------------    ERROR #1    -------------|")
    print("|       Time seriesC is not binary         |")
    print("|      (or the number of events=0)!        |")
    print("| Use CC.binarize() to preprocess seriesC  |")
    print("|------------------------------------------|")
    return()
  }
  if (tau < 0 || delT < 0) {
    print("|-------------    ERROR #2    -------------|")
    print("|    The offset (tau) or delta T (delT)    |")
    print("|              is negative.                |")
    print("|------------------------------------------|")
    return()
  }
  if (length(seriesA) != length(seriesB) || length(seriesA) != 
      length(seriesC) || length(seriesB) != length(seriesC)) {
    print("|-------------    ERROR #4    -------------|")
    print("|    lengths of series A, B and C differ   |")
    print("|------------------------------------------|")
    return()
  }
  if (!is.vector(seriesA) && !is.matrix(seriesA)) {
    print("|-------------    ERROR #5    -------------|")
    print("| series A is neither vector nor matrix    |")
    print("|------------------------------------------|")
    return()
  }
  if (!is.vector(seriesB) && !is.matrix(seriesB)) {
    print("|-------------    ERROR #5    -------------|")
    print("| series B is neither vector nor matrix    |")
    print("|------------------------------------------|")
    return()
  }
  if (!is.vector(seriesC) && !is.matrix(seriesC)) {
    print("|-------------    ERROR #5    -------------|")
    print("| series C is neither vector nor matrix    |")
    print("|------------------------------------------|")
    return()
  }
  if (any(is.na(seriesA)) && sigtest == "surrogate") {
    print("|-------------    ERROR #6    -------------|")
    print("|          seriesA contains NAs.           |")
    print("| surrogate test not allowed, use poisson  |")
    print("|------------------------------------------|")
    return()
  }
  if (any(is.na(seriesB)) && sigtest == "surrogate") {
    print("|-------------    ERROR #6    -------------|")
    print("|          seriesB contains NAs.           |")
    print("| surrogate test not allowed, use poisson  |")
    print("|------------------------------------------|")
    return()
  }
  if (any(is.na(seriesB)) && sigtest == "shuffle") {
    print("|-------------    ERROR #6    -------------|")
    print("|    seriesB or seriesA contains NAs.      |")
    print("| shuffle test not allowed, use poisson    |")
    print("|------------------------------------------|")
    return()
  }
  if (any(is.na(seriesC)) && sigtest == "surrogate") {
    print("|-------------    ERROR #6    -------------|")
    print("|          seriesC contains NAs.           |")
    print("| surrogate test not allowed, use poisson  |")
    print("|------------------------------------------|")
    return()
  }
  if (sum(is.na(seriesA)) == length(seriesA) || sum(is.na(seriesB)) == 
      length(seriesB) || sum(is.na(seriesC)) == length(seriesC)) {
    CA_out = list(NA, NA, NA, NA,NA, NA, NA, NA, NA,NA)
    names(CA_out) = c("NH precursor", "NH trigger", "p-value precursor", 
                      "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                      "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")
    return(CA_out)
  }
  bindata = rbind(seriesB, seriesC)
  rownames(bindata) = c("seriesB", "seriesC")
  Tlen = length(bindata[1, !is.na(bindata[1, ])])
  steplen = length(seriesB)
  seriesBC = matrix(0, 1, steplen)
  seriesBC[1, is.na(bindata[1, ])] = NA
  N_B = as.numeric(Tlen - table(bindata[1, ] == 1)[1])
  N_C = as.numeric(Tlen - table(bindata[2, ] == 1)[1])
  for (step in 1:steplen) {
    if (is.na(bindata[1, step])) {
      next
    }
    if (bindata[1, step] == 1) {
      start = step - tau.cond - (delT.cond)
      if (sym.cond == TRUE) {
        end = step - tau.cond + (delT.cond)
      }
      if (sym.cond == FALSE) {
        end = step - tau.cond
      }
      if (start < 1 & end >= 1) {
        start = 1
      }
      if (end < 1) {
        next
      }
      if (start > steplen) {
        next
      }
      if (end > steplen) {
        end = steplen
      }
      if (is.element(1, bindata[2, start:end])) {
        seriesBC[1, step] = 1
      }
    }
  }
  if (length(table(seriesBC)) < 2) {
    sig_testprec = TRUE
    sig_testtrigg = TRUE
    Pprec = 1
    Ptrigg = 1
    CRprec = 0
    CRtrigg = 0
    CA_out = list(sig_testprec, sig_testtrigg, Pprec, Ptrigg, 
                  CRprec, CRtrigg, 0, 0, NA,NA)
    print(paste0("No C-conditioned B-events found using delT.cond = ", 
                 delT.cond, "and tau.cond = ", tau.cond))
    names(CA_out) = c("NH precursor", "NH trigger", "p-value precursor", 
                      "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                      "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")
    return(CA_out)
  }
  seriesA[is.na(seriesBC)] = NA
  seriesBC[is.na(seriesA)] = NA
  bindata = rbind(seriesA, seriesBC)
  rownames(bindata) = c("seriesA", "seriesBC")
  Tlen = length(bindata[1, !is.na(bindata[1, ])])
  steplen = length(seriesA)
  N_A = as.numeric(Tlen - table(bindata[1, ] == 1)[1])
  N_B = as.numeric(Tlen - table(bindata[2, ] == 1)[1])
  K_prec = 0
  for (step in 1:steplen) {
    if (is.na(bindata[1, step])) {
      next
    }
    if (bindata[1, step] == 1) {
      start = step - tau - (delT)
      if (sym == TRUE) {
        end = step - tau + (delT)
      }
      if (sym == FALSE) {
        end = step - tau
      }
      if (start < 1 & end >= 1) {
        start = 1
      }
      if (end < 1) {
        next
      }
      if (start > steplen) {
        next
      }
      if (end > steplen) {
        end = steplen
      }
      if (is.element(1, bindata[2, start:end])) {
        K_prec = K_prec + 1
      }
    }
  }
  CRprec = K_prec/N_A
  K_trigg = 0
  for (step in 1:steplen) {
    if (is.na(bindata[2, step])) {
      next
    }
    if (bindata[2, step] == 1) {
      if (sym == TRUE) {
        start = step + tau - (delT)
      }
      if (sym == FALSE) {
        start = step + tau
      }
      end = step + tau + (delT)
      if (start < 1 & end >= 1) {
        start = 1
      }
      if (end < 1) {
        next
      }
      if (start > steplen) {
        next
      }
      if (end > steplen) {
        end = steplen
      }
      if (is.element(1, bindata[1, start:end])) {
        K_trigg = K_trigg + 1
      }
    }
  }
  CRtrigg = K_trigg/N_B
  if (sigtest == "poisson") {
    if (sym == FALSE) {
      Pprec = 0
      for (Ktmp in K_prec:N_A) {
        Ptmp = choose(N_A, Ktmp) * ((1 - ((1 - ((delT + 
                                                   1)/(Tlen - tau)))^N_B))^Ktmp) * (((1 - ((delT + 
                                                                                              1)/(Tlen - tau)))^N_B)^(N_A - Ktmp))
        Pprec = Pprec + Ptmp
        Eprec = 1 - (1 - delT/(Tlen - tau))^N_B
      }
      Ptrigg = 0
      for (Ktmp in K_trigg:N_B) {
        Ptmp = choose(N_B, Ktmp) * ((1 - ((1 - ((delT + 
                                                   1)/(Tlen - tau)))^N_A))^Ktmp) * (((1 - ((delT + 
                                                                                              1)/(Tlen - tau)))^N_A)^(N_B - Ktmp))
        Ptrigg = Ptrigg + Ptmp
        Etrigg = 1 - (1 - delT/(Tlen - tau))^N_A
      }
    }
    if (sym == TRUE) {
      Pprec = 0
      delTsym = delT + delT + 1
      for (Ktmp in K_prec:N_A) {
        Ptmp = choose(N_A, Ktmp) * ((1 - ((1 - (delTsym/(Tlen)))^N_B))^Ktmp) * 
          (((1 - (delTsym/(Tlen)))^N_B)^(N_A - Ktmp))
        Pprec = Pprec + Ptmp
        Eprec = 1 - (1 - delT/(Tlen - tau))^N_B
      }
      Ptrigg = 0
      for (Ktmp in K_trigg:N_B) {
        Ptmp = choose(N_B, Ktmp) * ((1 - ((1 - (delTsym/(Tlen)))^N_A))^Ktmp) * 
          (((1 - (delTsym/(Tlen)))^N_A)^(N_B - Ktmp))
        Ptrigg = Ptrigg + Ptmp
        Etrigg = 1 - (1 - delT/(Tlen - tau))^N_A
      }
    }
  }
  if (sigtest == "wt.surrogate") {
    if (N_A > 2 & N_B > 2){
    if (delT == 0) {
      delT = 1
    }
    seriesA = CC.ts2es(seriesA)
    seriesB = CC.ts2es(seriesB)
    span = (max(seriesA$span[1], seriesB$span[1]):min(seriesA$span[2], 
                                                      seriesB$span[2]))
    Tlen = length(span)
    seriesA_sort = sort(seriesA$es)
    wtA = rep(0, N_A)
    wtA[1] = seriesA_sort[1] - seriesA$span[1]
    for (i in 2:N_A) {
      wtA[i] = seriesA_sort[i] - seriesA_sort[i - 1]
    }
    seriesB_sort = sort(seriesB$es)
    wtB = rep(0, N_B)
    wtB[1] = seriesB_sort[1] - seriesB$span[1]
    for (i in 2:N_B) {
      wtB[i] = seriesB_sort[i] - seriesB_sort[i - 1]
    }
    surdist = matrix(NA, 2, reps)
    for (surno in 1:reps) {
      surA = sample(wtA, size = 1)
      for (i in 2:length(seriesA$es)) {
        tmp = surA[i - 1] + sample(wtA, size = 1)
        if (tmp > span[Tlen]) {
          break
        }
        else {
          surA[i] = tmp
        }
      }
      surB = sample(wtB, size = 1)
      for (i in 2:length(seriesB$es)) {
        tmp = surB[i - 1] + sample(wtB, size = 1)
        if (tmp > span[Tlen]) {
          break
        }
        else {
          surB[i] = tmp
        }
      }
      K_prec_sur = 0
      for (step_a in 1:N_A) {
        if (!is.element(surA[step_a], span)) {
          next
        }
        for (step_b in 1:N_B) {
          if (!is.element(surB[step_b], span)) {
            next
          }
          if (sym == FALSE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(0:delT))) {
              K_prec_sur = K_prec_sur + 1
              break
            }
          }
          if (sym == TRUE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(-delT:delT))) {
              K_prec_sur = K_prec_sur + 1
              break
            }
          }
        }
      }
      surdist[1, surno] = K_prec_sur/N_A
      K_trigg_sur = 0
      for (step_b in 1:N_B) {
        if (!is.element(surB[step_b], span)) {
          next
        }
        for (step_a in 1:N_A) {
          if (!is.element(surA[step_a], span)) {
            next
          }
          if (sym == FALSE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(0:delT))) {
              K_trigg_sur = K_trigg_sur + 1
              break
            }
          }
          if (sym == TRUE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(-delT:delT))) {
              K_trigg_sur = K_trigg_sur + 1
              break
            }
          }
        }
      }
      surdist[2, surno] = K_trigg_sur/N_B
    }
    Pprec = 1 - ecdf(surdist[1, ])(CRprec)
    Ptrigg = 1 - ecdf(surdist[2, ])(CRtrigg)
    Eprec = mean(surdist[1, ])
    Etrigg = mean(surdist[2, ])
    }else{
      Pprec = NaN
      Ptrigg = NaN
      Eprec = NaN
      Etrigg = NaN
      }
  }
  if (sigtest == "shuffle.surrogate") {
    if (N_A > 2 & N_B > 2){
    if (delT == 0) {
      delT = 1
    }
    surdist = matrix(NA, 2, reps)
    span = seq(1:steplen)
    for (surno in 1:reps) {
      surA = sample(span, size = N_A)
      surB = sample(span, size = N_B)
      K_prec_sur = 0
      for (step_a in 1:N_A) {
        for (step_b in 1:N_B) {
          if (sym == FALSE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(0:delT))) {
              K_prec_sur = K_prec_sur + 1
              break
            }
          }
          if (sym == TRUE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(-delT:delT))) {
              K_prec_sur = K_prec_sur + 1
              break
            }
          }
        }
      }
      surdist[1, surno] = K_prec_sur/N_A
      K_trigg_sur = 0
      for (step_b in 1:N_B) {
        for (step_a in 1:N_A) {
          if (sym == FALSE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(0:delT))) {
              K_trigg_sur = K_trigg_sur + 1
              break
            }
          }
          if (sym == TRUE) {
            if (is.element(((surA[step_a] - tau) - surB[step_b]), 
                           c(-delT:delT))) {
              K_trigg_sur = K_trigg_sur + 1
              break
            }
          }
        }
      }
      surdist[2, surno] = K_trigg_sur/N_B
    }
    Pprec = 1 - ecdf(surdist[1, ])(CRprec)
    Ptrigg = 1 - ecdf(surdist[2, ])(CRtrigg)
    Eprec = mean(surdist[1, ])
    Etrigg = mean(surdist[2, ])
    }else{
      Pprec = NaN
      Ptrigg = NaN
      Eprec = NaN
      Etrigg = NaN
      
    }
  }
  
  
  if(is.na(Pprec)){
    sig_testprec = NaN
    sig_testtrigg = NaN
  }else{
  sig_testprec = logical
  if (Pprec >= alpha) {
    sig_testprec = TRUE
  }
  else {
    sig_testprec = FALSE
  }
  sig_testtrigg = logical
  if (Ptrigg >= alpha) {
    sig_testtrigg = TRUE
  }
  else {
    sig_testtrigg = FALSE
  }
  }
  CA_out = list(sig_testprec, sig_testtrigg, Pprec, Ptrigg, 
                CRprec, CRtrigg, N_A, N_B, Eprec,Etrigg)
  names(CA_out) = c("NH precursor", "NH trigger", "p-value precursor", 
                    "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                    "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")
  return(CA_out)
}


dat_to_adj <- function(x){
  # x -> c1 - time; c2 - from (origin); c3 - to (destination)
  n <- max(x[,1])
  fls <- list()
  for (i in 1:n){
    fls[[i]] <- x[x[,1] == i  & x[,2] != x[,3],2:3]
    row.names(fls[[i]]) <- NULL
  }
  return(fls)
}

net_recip <- function(x, lag){
  # x - output from dat_to_adj(); lag - your time lag
  n <- length(x)
  out <- c()
  for (t in 1:(n-lag)){
    x1 <- x[[t]]
    x2 <- data.frame(fromDist = x[[t+lag]][,2],toDist = x[[t+lag]][,1])
    x1_bar <- nrow(x1)/74/73
    x2_bar <- nrow(x2)/74/73
    top_val <- 0
    bot_val1 <- 0
    bot_val2 <- 0
    for (i in 1:74){
      for (j in 1:74){
        if (i != j){
          a1 <- sum(x1[,1] == i & x1[,2] == j)
          a2 <- sum(x2[,1] == i & x2[,2] == j)
          top_val <- top_val + (a1 - x1_bar) * (a2 - x2_bar)
          bot_val1 <- bot_val1 + (a1 - x1_bar)^2
          bot_val2 <- bot_val2 + (a2 - x2_bar)^2
        }
      }
    }
    rho <- top_val / sqrt(bot_val1) / sqrt(bot_val2)
    out <- c(out, rho)
  }
  return(out)
}

uhuy <- function(idp_dat,conf_dat,rain_dat,tmax_dat,n,dt){
  no_idp <- colSums(idp_dat) > 0
  conf_idp_other_lag <- data.frame(matrix(0, ncol = n+1, nrow = 74))
  rain_idp_other_lag <- data.frame(matrix(0, ncol = n+1, nrow = 74))
  tmax_idp_other_lag <- data.frame(matrix(0, ncol = n+1, nrow = 74))
  
  rain_conf_other_lag <- data.frame(matrix(0, ncol = n+1, nrow = 74))
  tmax_conf_other_lag <- data.frame(matrix(0, ncol = n+1, nrow = 74))
  
  for (i in 1:74){
    for (t in 1:n){
      idp_dat2 <- idp_dat[,i]
      conf_dat2 <- conf_dat[,i]
      rain_dat2 <- rain_dat[,i]
      tmax_dat2 <- tmax_dat[,i]
      
      # conflict-migration
      if (sum(idp_dat2) * sum(conf_dat2) != 0){
        out_conf_idp <- eca_two(idp_dat2, conf_dat2,alpha=0.05,delT=dt,tau=t-1,
                                sigtest="poisson",reps = 1000)$`precursor coincidence rate`
      }else{out_conf_idp <- -1}
      # rainfall-migration
      if (sum(idp_dat2) * sum(rain_dat2) != 0){
        out_rain_idp <- eca_two(idp_dat2, rain_dat2,alpha=0.05,delT=dt,tau=t-1,
                                sigtest="poisson",reps = 1000)$`precursor coincidence rate`
      }else{out_rain_idp <- -1}
      # rainfall-conflict
      if (sum(conf_dat2) * sum(rain_dat2) != 0){
        out_rain_conf <- eca_two(conf_dat2, rain_dat2,alpha=0.05,delT=dt,tau=t-1,
                                 sigtest="poisson",reps = 1000)$`precursor coincidence rate`
      }else{out_rain_conf <- -1}
      # temp-migration
      if (sum(idp_dat2) * sum(tmax_dat2) != 0){
        out_tmax_idp <- eca_two(idp_dat2, tmax_dat2,alpha=0.05,delT=dt,tau=t-1,
                                sigtest="poisson",reps = 1000)$`precursor coincidence rate`
      }else{out_tmax_idp <- -1}
      # temp-conflict
      if (sum(conf_dat2) * sum(tmax_dat2) != 0){
        out_tmax_conf <- eca_two(conf_dat2, tmax_dat2,alpha=0.05,delT=dt,tau=t-1,
                                 sigtest="poisson",reps = 1000)$`precursor coincidence rate`
      }else{out_tmax_conf <- -1}
      
      conf_idp_other_lag[i,t] <- out_conf_idp
      rain_idp_other_lag[i,t] <- out_rain_idp
      rain_conf_other_lag[i,t] <- out_rain_conf
      tmax_idp_other_lag[i,t] <- out_tmax_idp
      tmax_conf_other_lag[i,t] <- out_tmax_conf
      
    }
    
    if (which.max(conf_idp_other_lag[i,1:n])-1 < 0) {
      conf_idp_other_lag[i,n+1] <- 0
    }else{conf_idp_other_lag[i,n+1] <- which.max(conf_idp_other_lag[i,1:n])-1}
    
    if (which.max(rain_idp_other_lag[i,1:n])-1 < 0) {
      rain_idp_other_lag[i,n+1] <- 0
    }else{rain_idp_other_lag[i,n+1] <- which.max(rain_idp_other_lag[i,1:n])-1}
    
    if (which.max(rain_conf_other_lag[i,1:n])-1 < 0) {
      rain_conf_other_lag[i,n+1] <- 0
    }else{rain_conf_other_lag[i,n+1] <- which.max(rain_conf_other_lag[i,1:n])-1}
    
    if (which.max(tmax_idp_other_lag[i,1:n])-1 < 0) {
      tmax_idp_other_lag[i,n+1] <- 0
    }else{tmax_idp_other_lag[i,n+1] <- which.max(tmax_idp_other_lag[i,1:n])-1}
    
    if (which.max(tmax_conf_other_lag[i,1:n])-1 < 0) {
      tmax_conf_other_lag[i,n+1] <- 0
    }else{tmax_conf_other_lag[i,n+1] <- which.max(tmax_conf_other_lag[i,1:n])-1}
    #
  }
  out <- list(conf_idp_other_lag,rain_idp_other_lag,rain_conf_other_lag,tmax_idp_other_lag,tmax_conf_other_lag)
  return(out)
}

real_eca <- function(idp_dat,conf_dat,rain_dat,tmax_dat,n,dt,lag_out,sigt){
  conf_idp_other <- data.frame(matrix(0, ncol = 10, nrow = 74))
  rain_idp_other <- data.frame(matrix(0, ncol = 10, nrow = 74))
  rain_conf_idp_other <- data.frame(matrix(0, ncol = 10, nrow = 74))
  tmax_idp_other <- data.frame(matrix(0, ncol = 10, nrow = 74))
  tmax_conf_idp_other <- data.frame(matrix(0, ncol = 10, nrow = 74))
  
  names(conf_idp_other) = c("NH precursor", "NH trigger", "p-value precursor", 
                            "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                            "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")
  names(rain_idp_other) = c("NH precursor", "NH trigger", "p-value precursor", 
                            "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                            "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")
  names(rain_conf_idp_other) = c("NH precursor", "NH trigger", "p-value precursor", 
                                 "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                                 "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")
  names(tmax_idp_other) = c("NH precursor", "NH trigger", "p-value precursor", 
                            "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                            "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")
  names(tmax_conf_idp_other) = c("NH precursor", "NH trigger", "p-value precursor", 
                                 "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                                 "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")

  conf_idp_other_lag <- lag_out[[1]]
  rain_idp_other_lag <- lag_out[[2]]
  rain_conf_other_lag <- lag_out[[3]]
  tmax_idp_other_lag <- lag_out[[4]]
  tmax_conf_other_lag <- lag_out[[5]]
  
  for (i in 1:74){
    # run ECA
    # climate-conflict
    if (sum(conf_dat[,i]) >= 2 & sum(rain_dat[,i]) >= 2){
      rain_conf <- eca_two(conf_dat[,i], rain_dat[,i],alpha=0.05,delT=dt,tau=rain_conf_other_lag[i,n+1],
                           sigtest=sigt,reps = 100)
    }else{rain_conf <- data.frame(t(c(-1,-1,-1,-1,-1,-1,1,1,-1,-1)))
    names(rain_conf) = c("NH precursor", "NH trigger", "p-value precursor", 
                         "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                         "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")}
    if (sum(conf_dat[,i]) >= 2 & sum(tmax_dat[,i]) >= 2){
      tmax_conf <- eca_two(conf_dat[,i], tmax_dat[,i],alpha=0.05,delT=dt,tau=tmax_conf_other_lag[i,n+1],
                           sigtest=sigt,reps = 100)
    }else{tmax_conf <- data.frame(t(c(-1,-1,-1,-1,-1,-1,1,1,-1,-1)))
    names(tmax_conf) = c("NH precursor", "NH trigger", "p-value precursor", 
                         "p-value trigger", "precursor coincidence rate", "trigger coincidence rate", 
                         "N precursor", "N trigger","Expected precursor rate","Expected trigger rate")}
    
    # conflict-migration
    if (sum(idp_dat[,i]) >= 2 & sum(conf_dat[,i])  >= 2){
      out_conf_idp <- eca_two(idp_dat[,i], conf_dat[,i],alpha=0.05,delT=dt,tau=conf_idp_other_lag[i,n+1],
                              sigtest=sigt,reps = 100)
    }else{out_conf_idp <- c(NaN,NaN,NaN,NaN,NaN,NaN,sum(idp_dat[,i]),sum(conf_dat[,i]),NaN,NaN)}
    # rainfall-migration
    if (sum(idp_dat[,i]) >= 2 & sum(rain_dat[,i]) >= 2){
      out_rain_idp <- eca_two(idp_dat[,i], rain_dat[,i],alpha=0.05,delT=dt,tau=rain_idp_other_lag[i,n+1],
                              sigtest=sigt,reps = 100)
    }else{out_rain_idp <- c(NaN,NaN,NaN,NaN,NaN,NaN,sum(idp_dat[,i]),sum(rain_dat[,i]),NaN,NaN)}
    # rainfall-conflict-migration
    if (sum(conf_dat[,i])>= 2 &sum(rain_dat[,i])>= 2 &sum(idp_dat[,i]) >= 2 & 
        rain_conf$`precursor coincidence rate` * rain_conf$`N precursor` >= 2 &
        rain_conf$`trigger coincidence rate` * rain_conf$`N trigger` >= 2){
      out_rain_conf_idp <- eca_three(idp_dat[,i], conf_dat[,i], rain_dat[,i],alpha=0.05,delT=dt,
                                     tau=conf_idp_other_lag[i,n+1],delT.cond=dt,tau.cond=rain_conf_other_lag[i,n+1],
                                     sigtest=sigt,reps = 100)
    }else{out_rain_conf_idp <- c(NaN,NaN,NaN,NaN,NaN,NaN,sum(idp_dat[,i]),sum(conf_dat[,i]),NaN,NaN)}
    # temp-migration
    if (sum(idp_dat[,i])>= 2 &sum(tmax_dat[,i]) >= 2){
      out_tmax_idp <- eca_two(idp_dat[,i], tmax_dat[,i],alpha=0.05,delT=dt,tau=tmax_idp_other_lag[i,n+1],
                              sigtest="wt.surrogate",reps = 100)
    }else{out_tmax_idp <- c(NaN,NaN,NaN,NaN,NaN,NaN,sum(idp_dat[,i]),sum(tmax_dat[,i]),NaN,NaN)}
    # temp-conflict-migration
    if (sum(conf_dat[,i])>= 2 &sum(tmax_dat[,i]) >= 2 &sum(idp_dat[,i]) >= 2  & 
        tmax_conf$`precursor coincidence rate` * tmax_conf$`N precursor` >= 2 &
        tmax_conf$`trigger coincidence rate` * tmax_conf$`N trigger` >= 2){
      out_tmax_conf_idp <- eca_three(idp_dat[,i], conf_dat[,i], tmax_dat[,i],alpha=0.05,delT=dt,
                                     tau=conf_idp_other_lag[i,n+1],delT.cond=dt,tau.cond=tmax_conf_other_lag[i,n+1],
                                     sigtest="wt.surrogate",reps = 100)
    }else{out_tmax_conf_idp <- c(NaN,NaN,NaN,NaN,NaN,NaN,sum(idp_dat[,i]),sum(conf_dat[,i]),NaN,NaN)}
    
    conf_idp_other[i,] <- out_conf_idp
    rain_idp_other[i,] <- out_rain_idp
    rain_conf_idp_other[i,] <- out_rain_conf_idp
    
    tmax_idp_other[i,] <- out_tmax_idp
    tmax_conf_idp_other[i,] <- out_tmax_conf_idp
  }
  
  out <- list(conf_idp_other,rain_idp_other,rain_conf_idp_other,tmax_idp_other,tmax_conf_idp_other)
  return(out)
}

eca_clean <- function(pois,wtsur,shuf){
  pois$`NH precursor`[is.nan(pois$`NH precursor`)] = -1
  wtsur$`NH precursor`[is.nan(wtsur$`NH precursor`)] = -1
  shuf$`NH precursor`[is.nan(shuf$`NH precursor`)] = -1
  
  stat <- ((pois$`NH precursor` == 0)| 
                      (wtsur$`NH precursor` == 0)| 
                      (shuf$`NH precursor` == 0))
  return(stat)
}

eca_plot <- function(shp,dt,type, occ = TRUE){
  ratefile_rain = paste("coinratemap_rain_",type,dt,".png",sep ="")
  occfile_rain = paste("coinnummap_rain_",type,dt,".png",sep ="")
  ratefile_tmax = paste("coinratemap_tmax_",type,dt,".png",sep ="")
  occfile_tmax = paste("coinnummap_tmax_",type,dt,".png",sep ="")
  
  ratename_rain <- c(paste("conf_idp_pre_",type,dt,sep =""),paste("rain_idp_pre_",type,dt,sep =""),
                     paste("rain_conf_idp_pre_",type,dt,sep =""))
  occname_rain <- c(paste("conf_idp_pre_n_",type,dt,sep =""),paste("rain_idp_pre_n_",type,dt,sep =""),
                    paste("rain_conf_idp_pre_n_",type,dt,sep =""))
  ratename_tmax <- c(paste("conf_idp_pre_",type,dt,sep =""),paste("tmax_idp_pre_",type,dt,sep =""),
                     paste("tmax_conf_idp_pre_",type,dt,sep =""))
  occname_tmax <- c(paste("conf_idp_pre_n_",type,dt,sep =""),paste("tmax_idp_pre_n_",type,dt,sep =""),
                    paste("tmax_conf_idp_pre_n_",type,dt,sep =""))
  
  png(ratefile_rain, width = 2000, height = 800, res = 300)
  plot1 <- spplot(shp, ratename_rain[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  plot2 <-  spplot(shp, ratename_rain[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   ))
  
  plot3 <-  spplot(shp, ratename_rain[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,ncol = 3)
  dev.off()
  
  png(ratefile_tmax, width = 2000, height = 800, res = 300)
  plot1 <- spplot(shp, ratename_tmax[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  
  plot2 <-  spplot(shp, ratename_tmax[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  plot3 <-  spplot(shp, ratename_tmax[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,ncol = 3)
  dev.off()
  if (occ == TRUE){
    png(occfile_rain, width = 2000, height = 800, res = 300)
    plot1 <- spplot(shp, occname_rain[1],
                    col.regions = rev(heat.colors(101)),
                    at = seq(-0.001,365,length=100),
                    #main=list(label="Conflict &",cex=1),
                    sp.layout = list(
                      list("sp.polygons",shp, first = TRUE, fill = "gray")
                    )) 
    plot2 <- spplot(shp, occname_rain[2],
                    col.regions = rev(heat.colors(101)),
                    at = seq(-0.001,365,length=100),
                    #main=list(label="rain-IDP prec., dT = 4",cex=1),
                    sp.layout = list(
                      list("sp.polygons",shp, first = TRUE, fill = "gray")
                    ))
    
    plot3 <- spplot(shp, occname_rain[3],
                    col.regions = rev(heat.colors(101)),
                    at = seq(-0.001,365,length=100),
                    #main=list(label="rain-conf-IDP prec., dT = 4",cex=1),
                    sp.layout = list(
                      list("sp.polygons",shp, first = TRUE, fill = "gray")
                    ))
    
    grid.arrange(plot1,plot2,plot3,ncol = 3)
    dev.off()
    
    png(occfile_tmax, width = 2000, height = 800, res = 300)
    plot1 <- spplot(shp, occname_tmax[1],
                    col.regions = rev(heat.colors(101)),
                    at = seq(-0.001,365,length=100),
                    #main=list(label="Conflict &",cex=1),
                    sp.layout = list(
                      list("sp.polygons",shp, first = TRUE, fill = "gray")
                    )) 
    
    plot2 <- spplot(shp, occname_tmax[2],
                    col.regions = rev(heat.colors(101)),
                    at = seq(-0.001,365,length=100),
                    #main=list(label="tmax-IDP prec., dT = 4",cex=1),
                    sp.layout = list(
                      list("sp.polygons",shp, first = TRUE, fill = "gray")
                    ))
    
    plot3 <- spplot(shp, occname_tmax[3],
                    col.regions = rev(heat.colors(101)),
                    at = seq(-0.001,365,length=100),
                    #main=list(label="tmax-conf-IDP prec., dT = 4",cex=1),
                    sp.layout = list(
                      list("sp.polygons",shp, first = TRUE, fill = "gray")
                    ))
    
    grid.arrange(plot1,plot2,plot3,ncol = 3)
    dev.off()
  }
}

corr_plot <- function(BMOther,eca_pois,dt,th,nan_ctr = FALSE,type){
  BMOther$B[BMOther$boundCon == 'fail'] = NA
  df.ci <- data.frame(b = BMOther$B, c = eca_pois[[1]]$`precursor coincidence rate`)
  df.di <- data.frame(b = BMOther$B, c = eca_pois[[2]]$`precursor coincidence rate`)
  df.dci <- data.frame(b = BMOther$B, c = eca_pois[[3]]$`precursor coincidence rate`)
  df.ti <- data.frame(b = BMOther$B, c = eca_pois[[4]]$`precursor coincidence rate`)
  df.tci <- data.frame(b = BMOther$B, c = eca_pois[[5]]$`precursor coincidence rate`)
  
  if (nan_ctr == TRUE) {  
    df.ci$c[is.nan(df.ci$c)] <- 0
    df.di$c[is.nan(df.di$c)] <- 0
    df.dci$c[is.nan(df.dci$c)] <- 0
    df.ti$c[is.nan(df.ti$c)] <- 0
    df.tci$c[is.nan(df.tci$c)] <- 0
  }

  tval.ci <- !is.nan(df.ci$c) & !is.na(df.ci$b)
  tval.di <- !is.nan(df.di$c) & !is.na(df.di$b)
  tval.dci <- !is.nan(df.dci$c) & !is.na(df.dci$b)
  tval.ti <- !is.nan(df.ti$c) & !is.na(df.ti$b)
  tval.tci <- !is.nan(df.tci$c) & !is.na(df.tci$b)
  
  df.ci <- df.ci[tval.ci == 1,]
  df.di <- df.di[tval.di == 1,]
  df.dci <- df.dci[tval.dci == 1,]
  df.ti <- df.ti[tval.ti == 1,]
  df.tci <- df.tci[tval.tci == 1,]
  
  if (cor.test(df.ci$b,df.ci$c)$p.value < th){
    ggscatter(df.ci, x = "b", y = "c", size = 1,
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Burstiness", ylab = "Coincidence rate",
              add.params = list(color = "red"),xlim = c(-1,1),ylim = c(0,1),
              cor.coeff.args = list(label.x = 0.45,label.y = 0.75, label.sep = "\n"))
  }else{
    ggscatter(df.ci, x = "b", y = "c", size = 1,
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Burstiness", ylab = "Coincidence rate",
              add.params = list(color = "gray"),xlim = c(-1,1),ylim = c(0,1),
              cor.coeff.args = list(label.x = 0.45,label.y = 0.75, label.sep = "\n"))
  }
  ggsave(paste("rain_",type,"_ci",dt,".png",sep = ""), height = 2, width = 3.6)
  
  if (cor.test(df.di$b,df.di$c)$p.value < th){
    ggscatter(df.di, x = "b", y = "c", size = 1,
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Burstiness", ylab = "Coincidence rate",
              add.params = list(color = "gray"),xlim = c(-1,1),ylim = c(0,1),
              cor.coeff.args = list(label.x = 0.45,label.y = 0.75, label.sep = "\n"))
  }else{
    ggscatter(df.di, x = "b", y = "c", size = 1,
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Burstiness", ylab = "Coincidence rate",
              add.params = list(color = "gray"),xlim = c(-1,1),ylim = c(0,1),
              cor.coeff.args = list(label.x = 0.45,label.y = 0.75, label.sep = "\n"))
  }
  ggsave(paste("rain_",type,"_di",dt,".png",sep = ""), height = 2, width = 3.6)
  
  if (cor.test(df.dci$b,df.dci$c)$p.value < th){
    ggscatter(df.dci, x = "b", y = "c", size = 1,
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Burstiness", ylab = "Coincidence rate",
              add.params = list(color = "red"),xlim = c(-1,1),ylim = c(0,1),
              cor.coeff.args = list(label.x = 0.45,label.y = 0.75, label.sep = "\n"))
  }else{
    ggscatter(df.dci, x = "b", y = "c", size = 1,
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Burstiness", ylab = "Coincidence rate",
              add.params = list(color = "gray"),xlim = c(-1,1),ylim = c(0,1),
              cor.coeff.args = list(label.x = 0.45,label.y = 0.75, label.sep = "\n"))
  }
  ggsave(paste("rain_",type,"_dci",dt,".png",sep = ""), height = 2, width = 3.6)
  
  
  if (cor.test(df.ti$b,df.ti$c)$p.value < th){
    ggscatter(df.ti, x = "b", y = "c", size = 1,
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Burstiness", ylab = "Coincidence rate",
              add.params = list(color = "red"),xlim = c(-1,1),ylim = c(0,1),
              cor.coeff.args = list(label.x = 0.45,label.y = 0.75, label.sep = "\n"))
  }else{
    ggscatter(df.ti, x = "b", y = "c", size = 1,
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Burstiness", ylab = "Coincidence rate",
              add.params = list(color = "gray"),xlim = c(-1,1),ylim = c(0,1),
              cor.coeff.args = list(label.x = 0.45,label.y = 0.75, label.sep = "\n"))
  }
  ggsave(paste("tmax_",type,"_ti",dt,".png",sep = ""), height = 2, width = 3.6)
  
  if (cor.test(df.tci$b,df.tci$c)$p.value < th){
    ggscatter(df.tci, x = "b", y = "c", size = 1,
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Burstiness", ylab = "Coincidence rate",
              add.params = list(color = "red"),xlim = c(-1,1),ylim = c(0,1),
              cor.coeff.args = list(label.x = 0.45,label.y = 0.75, label.sep = "\n"))
  }else{
    ggscatter(df.tci, x = "b", y = "c", size = 1,
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Burstiness", ylab = "Coincidence rate",
              add.params = list(color = "gray"),xlim = c(-1,1),ylim = c(0,1),
              cor.coeff.args = list(label.x = 0.45,label.y = 0.75, label.sep = "\n"))
  }
  ggsave(paste("tmax_",type,"_tci",dt,".png",sep = ""), height = 2, width = 3.6)
}

agg_cal <- function(shp,timeDistOther){
  conf_idp_r <- shp$conf_idp_pre_other1
  conf_idp_r[is.nan(conf_idp_r)] <- 0
  agg_conf_idp <- sum(conf_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_conf_idp2 <- sum(conf_idp_r * colSums(timeDistOther)) / 74
  
  rain_idp_r <- shp$rain_idp_pre_other1
  rain_idp_r[is.nan(rain_idp_r)] <- 0
  agg_rain_idp <- sum(rain_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_rain_idp2 <- sum(rain_idp_r * colSums(timeDistOther)) / 74
  
  rain_conf_idp_r <- shp$rain_conf_idp_pre_other1
  rain_conf_idp_r[is.nan(rain_conf_idp_r)] <- 0
  agg_rain_conf_idp <- sum(rain_conf_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_rain_conf_idp2 <- sum(rain_conf_idp_r * colSums(timeDistOther)) / 74
  
  tmax_idp_r <- shp$tmax_idp_pre_other1
  tmax_idp_r[is.nan(tmax_idp_r)] <- 0
  agg_tmax_idp <- sum(tmax_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_tmax_idp2 <- sum(tmax_idp_r * colSums(timeDistOther)) / 74
  
  tmax_conf_idp_r <- shp$tmax_conf_idp_pre_other1
  tmax_conf_idp_r[is.nan(tmax_conf_idp_r)] <- 0
  agg_tmax_conf_idp <- sum(tmax_conf_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_tmax_conf_idp2 <- sum(tmax_conf_idp_r * colSums(timeDistOther)) / 74
  
  out_r1 <- c(agg_conf_idp,agg_rain_idp,agg_rain_conf_idp,agg_tmax_idp,agg_tmax_conf_idp)
  out_n1 <- c(agg_conf_idp2,agg_rain_idp2,agg_rain_conf_idp2,agg_tmax_idp2,agg_tmax_conf_idp2)
  
  conf_idp_r <- shp$conf_idp_pre_other4
  conf_idp_r[is.nan(conf_idp_r)] <- 0
  agg_conf_idp <- sum(conf_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_conf_idp2 <- sum(conf_idp_r * colSums(timeDistOther)) / 74
  
  rain_idp_r <- shp$rain_idp_pre_other4
  rain_idp_r[is.nan(rain_idp_r)] <- 0
  agg_rain_idp <- sum(rain_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_rain_idp2 <- sum(rain_idp_r * colSums(timeDistOther)) / 74
  
  rain_conf_idp_r <- shp$rain_conf_idp_pre_other4
  rain_conf_idp_r[is.nan(rain_conf_idp_r)] <- 0
  agg_rain_conf_idp <- sum(rain_conf_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_rain_conf_idp2 <- sum(rain_conf_idp_r * colSums(timeDistOther)) / 74
  
  tmax_idp_r <- shp$tmax_idp_pre_other4
  tmax_idp_r[is.nan(tmax_idp_r)] <- 0
  agg_tmax_idp <- sum(tmax_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_tmax_idp2 <- sum(tmax_idp_r * colSums(timeDistOther)) / 74
  
  tmax_conf_idp_r <- shp$tmax_conf_idp_pre_other4
  tmax_conf_idp_r[is.nan(tmax_conf_idp_r)] <- 0
  agg_tmax_conf_idp <- sum(tmax_conf_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_tmax_conf_idp2 <- sum(tmax_conf_idp_r * colSums(timeDistOther)) / 74
  
  out_r4 <- c(agg_conf_idp,agg_rain_idp,agg_rain_conf_idp,agg_tmax_idp,agg_tmax_conf_idp)
  out_n4 <- c(agg_conf_idp2,agg_rain_idp2,agg_rain_conf_idp2,agg_tmax_idp2,agg_tmax_conf_idp2)
  
  conf_idp_r <- shp$conf_idp_pre_other13
  conf_idp_r[is.nan(conf_idp_r)] <- 0
  agg_conf_idp <- sum(conf_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_conf_idp2 <- sum(conf_idp_r * colSums(timeDistOther)) / 74
  
  rain_idp_r <- shp$rain_idp_pre_other13
  rain_idp_r[is.nan(rain_idp_r)] <- 0
  agg_rain_idp <- sum(rain_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_rain_idp2 <- sum(rain_idp_r * colSums(timeDistOther)) / 74
  
  rain_conf_idp_r <- shp$rain_conf_idp_pre_other13
  rain_conf_idp_r[is.nan(rain_conf_idp_r)] <- 0
  agg_rain_conf_idp <- sum(rain_conf_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_rain_conf_idp2 <- sum(rain_conf_idp_r * colSums(timeDistOther)) / 74
  
  tmax_idp_r <- shp$tmax_idp_pre_other13
  tmax_idp_r[is.nan(tmax_idp_r)] <- 0
  agg_tmax_idp <- sum(tmax_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_tmax_idp2 <- sum(tmax_idp_r * colSums(timeDistOther)) / 74
  
  tmax_conf_idp_r <- shp$tmax_conf_idp_pre_other13
  tmax_conf_idp_r[is.nan(tmax_conf_idp_r)] <- 0
  agg_tmax_conf_idp <- sum(tmax_conf_idp_r * colSums(timeDistOther)) / sum(colSums(timeDistOther))
  agg_tmax_conf_idp2 <- sum(tmax_conf_idp_r * colSums(timeDistOther)) / 74
  
  out_r13 <- c(agg_conf_idp,agg_rain_idp,agg_rain_conf_idp,agg_tmax_idp,agg_tmax_conf_idp)
  out_n13 <- c(agg_conf_idp2,agg_rain_idp2,agg_rain_conf_idp2,agg_tmax_idp2,agg_tmax_conf_idp2)
  
  print(rbind(out_r1,out_r4,out_r13))
  print(rbind(out_n1,out_n4,out_n13))
  
}


si_exp_plot <- function(shp,dt){
  expfile_rain_pois = paste("expmap_rain_pois",dt,".png",sep ="")
  expfile_tmax_pois = paste("expmap_tmax_pois",dt,".png",sep ="")
  expfile_rain_wtsur = paste("expmap_rain_wtsur",dt,".png",sep ="")
  expfile_tmax_wtsur = paste("expmap_tmax_wtsur",dt,".png",sep ="")
  expfile_rain_shuf = paste("expmap_rain_shuf",dt,".png",sep ="")
  expfile_tmax_shuf = paste("expmap_tmax_shuf",dt,".png",sep ="")
  
  expname_rain_pois <- c(paste("expected_conf_idp_exp_pois",dt,sep =""),
                         paste("expected_rain_idp_exp_pois",dt,sep =""),
                         paste("expected_rain_conf_idp_exp_pois",dt,sep =""))
  expname_tmax_pois <- c(paste("expected_conf_idp_exp_pois",dt,sep =""),
                         paste("expected_tmax_idp_exp_pois",dt,sep =""),
                    paste("expected_tmax_conf_idp_exp_pois",dt,sep =""))
  expname_rain_wtsur <- c(paste("expected_conf_idp_exp_wtsur",dt,sep =""),
                         paste("expected_rain_idp_exp_wtsur",dt,sep =""),
                         paste("expected_rain_conf_idp_exp_wtsur",dt,sep =""))
  expname_tmax_wtsur <- c(paste("expected_conf_idp_exp_wtsur",dt,sep =""),
                         paste("expected_tmax_idp_exp_wtsur",dt,sep =""),
                         paste("expected_tmax_conf_idp_exp_wtsur",dt,sep =""))
  expname_rain_shuf <- c(paste("expected_conf_idp_exp_shuf",dt,sep =""),
                         paste("expected_rain_idp_exp_shuf",dt,sep =""),
                         paste("expected_rain_conf_idp_exp_shuf",dt,sep =""))
  expname_tmax_shuf <- c(paste("expected_conf_idp_exp_shuf",dt,sep =""),
                         paste("expected_tmax_idp_exp_shuf",dt,sep =""),
                         paste("expected_tmax_conf_idp_exp_shuf",dt,sep =""))
  
  # 1
  png(expfile_rain_pois, width = 2000, height = 800, res = 300)
  plot1 <- spplot(shp, expname_rain_pois[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  plot2 <-  spplot(shp, expname_rain_pois[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   ))
  
  plot3 <-  spplot(shp, expname_rain_pois[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,ncol = 3)
  dev.off()
  
  # 2
  png(expfile_tmax_pois, width = 2000, height = 800, res = 300)
  plot1 <- spplot(shp, expname_tmax_pois[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  plot2 <-  spplot(shp, expname_tmax_pois[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   ))
  
  plot3 <-  spplot(shp, expname_tmax_pois[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,ncol = 3)
  dev.off()
  
  #3
  png(expfile_rain_wtsur, width = 2000, height = 800, res = 300)
  plot1 <- spplot(shp, expname_rain_wtsur[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  plot2 <-  spplot(shp, expname_rain_wtsur[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   ))
  
  plot3 <-  spplot(shp, expname_rain_wtsur[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,ncol = 3)
  dev.off()
  
  # 4
  png(expfile_tmax_wtsur, width = 2000, height = 800, res = 300)
  plot1 <- spplot(shp, expname_tmax_wtsur[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  plot2 <-  spplot(shp, expname_tmax_wtsur[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   ))
  
  plot3 <-  spplot(shp, expname_tmax_wtsur[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,ncol = 3)
  dev.off()
  
  # 5
  png(expfile_rain_shuf, width = 2000, height = 800, res = 300)
  plot1 <- spplot(shp, expname_rain_shuf[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  plot2 <-  spplot(shp, expname_rain_shuf[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   ))
  
  plot3 <-  spplot(shp, expname_rain_shuf[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,ncol = 3)
  dev.off()
  
  # 6
  png(expfile_tmax_shuf, width = 2000, height = 800, res = 300)
  plot1 <- spplot(shp, expname_tmax_shuf[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  plot2 <-  spplot(shp, expname_tmax_shuf[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   ))
  
  plot3 <-  spplot(shp, expname_tmax_shuf[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,ncol = 3)
  dev.off()
  }

si_pval_plot <- function(shp,dt){
  pvalfile_pois = paste("pvalmap_pois",dt,".png",sep ="")
  pvalfile_wtsur = paste("pvalmap_wtsur",dt,".png",sep ="")
  pvalfile_shuf = paste("pvalmap_shuf",dt,".png",sep ="")
  pvalfile = paste("pvalmap",dt,".png",sep ="")

  pvalname_pois <- c(paste("conf_idp_pre_pv_pois",dt,sep =""),
                         paste("rain_idp_pre_pv_pois",dt,sep =""),
                         paste("rain_conf_idp_pre_pv_pois",dt,sep =""),
                     paste("tmax_idp_pre_pv_pois",dt,sep =""),
                     paste("tmax_conf_idp_pre_pv_pois",dt,sep =""))
  pvalname_wtsur <- c(paste("conf_idp_pre_pv_wtsur",dt,sep =""),
                         paste("rain_idp_pre_pv_wtsur",dt,sep =""),
                         paste("rain_conf_idp_pre_pv_wtsur",dt,sep =""),
                      paste("tmax_idp_pre_pv_wtsur",dt,sep =""),
                      paste("tmax_conf_idp_pre_pv_wtsur",dt,sep =""))

  pvalname_shuf <- c(paste("conf_idp_pre_pv_shuf",dt,sep =""),
                          paste("rain_idp_pre_pv_shuf",dt,sep =""),
                          paste("rain_conf_idp_pre_pv_shuf",dt,sep =""),
                     paste("tmax_idp_pre_pv_shuf",dt,sep =""),
                     paste("tmax_conf_idp_pre_pv_shuf",dt,sep =""))
  pvalname <- c(paste("conf_idp_pre_pv_other",dt,sep =""),
                          paste("rain_idp_pre_pv_other",dt,sep =""),
                          paste("rain_conf_idp_pre_pv_other",dt,sep =""),
                paste("tmax_idp_pre_pv_other",dt,sep =""),
                paste("tmax_conf_idp_pre_pv_other",dt,sep =""))

  # 1
  png(pvalfile_pois, width = 3500, height = 800, res = 300)
  plot1 <- spplot(shp, pvalname_pois[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  plot2 <-  spplot(shp, pvalname_pois[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   ))
  
  plot3 <-  spplot(shp, pvalname_pois[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  plot4 <-  spplot(shp, pvalname_pois[4],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  plot5 <-  spplot(shp, pvalname_pois[5],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,plot4,plot5,ncol = 5)
  dev.off()
  
  # 3
  png(pvalfile_wtsur, width = 3500, height = 800, res = 300)
  plot1 <- spplot(shp, pvalname_wtsur[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  plot2 <-  spplot(shp, pvalname_wtsur[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   ))
  
  plot3 <-  spplot(shp, pvalname_wtsur[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  plot4 <-  spplot(shp, pvalname_wtsur[4],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  plot5 <-  spplot(shp, pvalname_wtsur[5],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,plot4,plot5,ncol = 5)
  dev.off()
  
  # 5
  png(pvalfile_shuf, width = 3500, height = 800, res = 300)
  plot1 <- spplot(shp, pvalname_shuf[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  #main=list(label="Conflict &",cex=1),
                  sp.layout = list(
                    list("sp.polygons",shp, first = TRUE, fill = "gray")
                  )) 
  
  plot2 <-  spplot(shp, pvalname_shuf[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   ))
  
  plot3 <-  spplot(shp, pvalname_shuf[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  plot4 <-  spplot(shp, pvalname_shuf[4],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  plot5 <-  spplot(shp, pvalname_shuf[5],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   #main=list(label="Conflict &",cex=1),
                   sp.layout = list(
                     list("sp.polygons",shp, first = TRUE, fill = "gray")
                   )) 
  
  grid.arrange(plot1,plot2,plot3,plot4,plot5,ncol = 5)
  dev.off()
  
  # 7
  png(pvalfile, width = 3500, height = 800, res = 300)
  plot1 <- spplot(shp, pvalname[1],
                  col.regions = rev(heat.colors(101)),
                  at = seq(-0.001,1.001,length=100),
                  lwd = 1,
                  colorkey=FALSE) 
  
  plot2 <-  spplot(shp, pvalname[2],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   colorkey=FALSE)
  
  plot3 <-  spplot(shp, pvalname[3],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   colorkey=FALSE) 
  plot4 <-  spplot(shp, pvalname[4],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   colorkey=FALSE) 
  plot5 <-  spplot(shp, pvalname[5],
                   col.regions = rev(heat.colors(101)),
                   at = seq(-0.001,1.001,length=100),
                   lwd = 1,
                   colorkey=FALSE) 
  
  grid.arrange(plot1,plot2,plot3,plot4,plot5,ncol = 5)
  dev.off()
  
}

lag_eca <- function(idp_dat,conf_dat,rain_dat,tmax_dat,n,dt){
  no_idp <- colSums(idp_dat) > 0
  conf_idp_other_lag <- data.frame(matrix(0, ncol = n+1, nrow = 74))
  rain_idp_other_lag <- data.frame(matrix(0, ncol = n+1, nrow = 74))
  tmax_idp_other_lag <- data.frame(matrix(0, ncol = n+1, nrow = 74))
  
  rain_conf_other_lag <- data.frame(matrix(0, ncol = n+1, nrow = 74))
  tmax_conf_other_lag <- data.frame(matrix(0, ncol = n+1, nrow = 74))
  
  for (i in 1:74){
    for (t in 1:n){
      idp_dat2 <- idp_dat[,i]
      conf_dat2 <- conf_dat[,i]
      rain_dat2 <- rain_dat[,i]
      tmax_dat2 <- tmax_dat[,i]
      
      # conflict-migration
      if (sum(idp_dat2) * sum(conf_dat2) != 0){
        out_conf_idp <- eca_two(idp_dat2, conf_dat2,alpha=0.05,delT=dt,tau=t-1,
                                sigtest="poisson",reps = 1000)$`precursor coincidence rate`
      }else{out_conf_idp <- -1}
      # rainfall-migration
      if (sum(idp_dat2) * sum(rain_dat2) != 0){
        out_rain_idp <- eca_two(idp_dat2, rain_dat2,alpha=0.05,delT=dt,tau=t-1,
                                sigtest="poisson",reps = 1000)$`precursor coincidence rate`
      }else{out_rain_idp <- -1}
      # rainfall-conflict
      if (sum(conf_dat2) * sum(rain_dat2) != 0){
        out_rain_conf <- eca_two(conf_dat2, rain_dat2,alpha=0.05,delT=dt,tau=t-1,
                                 sigtest="poisson",reps = 1000)$`precursor coincidence rate`
      }else{out_rain_conf <- -1}
      # temp-migration
      if (sum(idp_dat2) * sum(tmax_dat2) != 0){
        out_tmax_idp <- eca_two(idp_dat2, tmax_dat2,alpha=0.05,delT=dt,tau=t-1,
                                sigtest="poisson",reps = 1000)$`precursor coincidence rate`
      }else{out_tmax_idp <- -1}
      # temp-conflict
      if (sum(conf_dat2) * sum(tmax_dat2) != 0){
        out_tmax_conf <- eca_two(conf_dat2, tmax_dat2,alpha=0.05,delT=dt,tau=t-1,
                                 sigtest="poisson",reps = 1000)$`precursor coincidence rate`
      }else{out_tmax_conf <- -1}
      
      conf_idp_other_lag[i,t] <- out_conf_idp
      rain_idp_other_lag[i,t] <- out_rain_idp
      rain_conf_other_lag[i,t] <- out_rain_conf
      tmax_idp_other_lag[i,t] <- out_tmax_idp
      tmax_conf_other_lag[i,t] <- out_tmax_conf
      
    }
    
    if (which.max(conf_idp_other_lag[i,1:n])-1 < 0) {
      conf_idp_other_lag[i,n+1] <- 0
    }else{conf_idp_other_lag[i,n+1] <- which.max(conf_idp_other_lag[i,1:n])-1}
    
    if (which.max(rain_idp_other_lag[i,1:n])-1 < 0) {
      rain_idp_other_lag[i,n+1] <- 0
    }else{rain_idp_other_lag[i,n+1] <- which.max(rain_idp_other_lag[i,1:n])-1}
    
    if (which.max(rain_conf_other_lag[i,1:n])-1 < 0) {
      rain_conf_other_lag[i,n+1] <- 0
    }else{rain_conf_other_lag[i,n+1] <- which.max(rain_conf_other_lag[i,1:n])-1}
    
    if (which.max(tmax_idp_other_lag[i,1:n])-1 < 0) {
      tmax_idp_other_lag[i,n+1] <- 0
    }else{tmax_idp_other_lag[i,n+1] <- which.max(tmax_idp_other_lag[i,1:n])-1}
    
    if (which.max(tmax_conf_other_lag[i,1:n])-1 < 0) {
      tmax_conf_other_lag[i,n+1] <- 0
    }else{tmax_conf_other_lag[i,n+1] <- which.max(tmax_conf_other_lag[i,1:n])-1}
    #
  }
  out <- list(conf_idp_other_lag,rain_idp_other_lag,rain_conf_other_lag,tmax_idp_other_lag,tmax_conf_other_lag)
  return(out)
}
