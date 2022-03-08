GetArgmaxId=function(calib_object=NULL, S=NULL){
  if ((is.null(calib_object))&(is.null(S))){
    stop("Please provide calib_object or S.")
  }
  
  if (is.null(S)){
    argmax_id=matrix(NA, nrow=ncol(calib_object$Lambda), ncol=2)
    for (block_id in 1:ncol(calib_object$Lambda)){
      if (ncol(calib_object$Lambda)==1){
        myS=calib_object$S
      } else {
        myS=calib_object$S_blocks[,block_id,drop=FALSE]
      }
      myS[is.na(myS)]=0
      myid=which.max(myS[,1])
      argmax_id[block_id,]=c(myid, which(calib_object$params$pi_list==calib_object$P[myid,block_id]))
    }
  } else {
    argmax_id=matrix(NA, nrow=1, ncol=2)
    myS=apply(S,1,max,na.rm=TRUE)
    myS[is.na(myS)]=0
    myid=which.max(myS)
    argmax_id[1,]=c(myid, max(which(S[myid,]==myS[myid])))
  }
  colnames(argmax_id)=c("lambda_id","pi_id")
  return(argmax_id)
}


GetArgmax=function(calib_object){
  argmax=matrix(NA, nrow=ncol(calib_object$Lambda), ncol=2)
  for (block_id in 1:ncol(calib_object$Lambda)){
    if (ncol(calib_object$Lambda)==1){
      myS=calib_object$S
    } else {
      myS=calib_object$S_blocks[,block_id,drop=FALSE]
    }
    myS[is.na(myS)]=0
    myid=which.max(myS[,1])
    argmax[block_id,]=c(calib_object$Lambda[myid,block_id], calib_object$P[myid,block_id])
  }
  colnames(argmax)=c("lambda","pi")
  return(argmax)
}


GetLambdaPath=function(lmax, lmin, cardinal=100){
  return(seq(sqrt(lmax),sqrt(lmin),length.out=cardinal)^2)
}


MakeBlockLambdaGrid=function(Lambda, lambda_other_blocks=NULL){
  if ((is.null(lambda_other_blocks))&(!is.vector(Lambda))){
    Lambda_blocks=Lambda
    Sequential_template=matrix(TRUE, ncol=ncol(Lambda), nrow=nrow(Lambda))
  } else {
    # Create Lambda grid matrix with nblocks columns
    if (!is.null(lambda_other_blocks)){
      nblocks=length(lambda_other_blocks)
    } else {
      lambda_other_blocks=1
      nblocks=1
    }
    Lambda_blocks=NULL
    if (is.vector(Lambda)){
      Sequential_template=matrix(FALSE, nrow=nblocks*length(Lambda), ncol=nblocks)
    } else {
      Sequential_template=matrix(FALSE, nrow=nblocks*nrow(Lambda), ncol=nblocks)
    }
    for (block_id in 1:nblocks){
      if (!is.vector(Lambda)){
        tmpLambda=Lambda[,block_id]
      } else {
        tmpLambda=Lambda
      }
      Lambda_blocks=cbind(Lambda_blocks, rep(lambda_other_blocks[block_id], nblocks*length(tmpLambda)))
      Lambda_blocks[(length(tmpLambda)*(block_id-1)+1):(length(tmpLambda)*(block_id)), block_id]=tmpLambda
      Sequential_template[(length(tmpLambda)*(block_id-1)+1):(length(tmpLambda)*(block_id)), block_id]=TRUE
    }
  }
  
  # # Remove duplicated rows 
  # rmids=which(!duplicated(apply(Lambda_blocks, 1, paste, collapse="-")))
  # Lambda_blocks=Lambda_blocks[rmids,,drop=FALSE]
  # Sequential_template=Sequential_template[rmids,,drop=FALSE]
  
  return(list(Lambda=Lambda_blocks, Sequential_template=Sequential_template))
}


GetContrast=function(mat, digits=3){
  return(length(unique(round(as.vector(abs(mat)), digits=digits))))
}


GetSubsample=function(ydata, family="gaussian", tau=0.5, resampling_method="subsampling", resampling_function=NULL){
  # resampling_function has to be a function of tau
  # argument resampling_method is ignored if resampling_function is provided
  if (resampling_method=="subsampling"){
    if (family=="gaussian"){
      if (is.null(resampling_function)){
        s=sample(length(ydata), size=tau*length(ydata))
      } else {
        s=resampling_function(tau=tau)
      }
    }
    if (family=="mgaussian"){
      if (is.null(resampling_function)){
        s=sample(nrow(ydata), size=tau*nrow(ydata))
      } else {
        s=resampling_function(tau=tau)
      }
    }
    if (family=="binomial"){
      if (is.null(resampling_function)){
        s0=sample(which(ydata=="0"), size=tau*sum(ydata=="0"))
        s1=sample(which(ydata=="1"), size=tau*sum(ydata=="1"))
        s=c(s0,s1)
      } else {
        s=resampling_function(tau=tau)
      }
    }
    if (family=="cox"){
      if (is.null(resampling_function)){
        s0=sample(which(ydata[,2]=="0"), size=tau*sum(ydata[,2]=="0"))
        s1=sample(which(ydata[,2]=="1"), size=tau*sum(ydata[,2]=="1"))
        s=c(s0,s1)
      } else {
        s=resampling_function(tau=tau)
      }
    }
  }
  if (resampling_method=="bootstrap"){
    if (family=="gaussian"){
      if (is.null(resampling_function)){
        s=sample(length(ydata), size=length(ydata), replace=TRUE)
      } else {
        s=resampling_function(tau=tau)
      }
    }
    if (family=="binomial"){
      if (is.null(resampling_function)){
        s0=sample(which(ydata=="0"), size=sum(ydata=="0"), replace=TRUE)
        s1=sample(which(ydata=="1"), size=sum(ydata=="1"), replace=TRUE)
        s=c(s0,s1)
      } else {
        s=resampling_function(tau=tau)
      }
    }
    if (family=="cox"){
      if (is.null(resampling_function)){
        s0=sample(which(ydata[,2]=="0"), size=sum(ydata[,2]=="0"), replace=TRUE)
        s1=sample(which(ydata[,2]=="1"), size=sum(ydata[,2]=="1"), replace=TRUE)
        s=c(s0,s1)
      } else {
        s=resampling_function(tau=tau)
      }
    }
  }
  return(s)
}


GetColorBar=function(colours, yspacing=1, yvalues=NULL, ylab="", mar=c(0,7,0,0), 
                     filename=NULL, width=3, height=8){
  if (is.null(yvalues)){
    yvalues=seq(1,length(colours))
  }
  if (!is.null(filename)){
    pdf(filename, width=width, height=height)
  }
  par(mar=mar)
  plot(NA,xlim=c(0,1), ylim=c(0,length(colours)), xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
  for (k in 1:length(colours)){
    polygon(x=c(0,1,1,0), y=c(k-1,k-1,k,k), col=colours[k], border=NA)
  }
  xseq=seq(1,length(yvalues),by=yspacing)
  axis(side=2, at=xseq, labels=yvalues[xseq], las=1, cex.axis=2)
  mtext(text=ylab, side=2, line=4, cex=4)
  if (!is.null(filename)){
    dev.off()
  }
}


GetBlockMatrix=function(pk){
  nblocks=sum(upper.tri(matrix(NA, ncol=length(pk), nrow=length(pk)), diag=TRUE))
  blocks=matrix(NA, nrow=length(pk), ncol=length(pk))
  blocks[upper.tri(blocks, diag=TRUE)]=1:nblocks
  
  mybreaks=c(0,cumsum(pk))
  bigblocks=matrix(ncol=sum(pk), nrow=sum(pk))
  row_id_start=matrix(mybreaks[row(blocks)],ncol=length(pk))+1
  row_id_end=matrix(mybreaks[row(blocks)+1],ncol=length(pk))
  col_id_start=matrix(mybreaks[col(blocks)],ncol=length(pk))+1
  col_id_end=matrix(mybreaks[col(blocks)+1],ncol=length(pk))
  
  row_id_start=row_id_start[upper.tri(row_id_start, diag=TRUE)]
  row_id_end=row_id_end[upper.tri(row_id_end, diag=TRUE)]
  col_id_start=col_id_start[upper.tri(col_id_start, diag=TRUE)]
  col_id_end=col_id_end[upper.tri(col_id_end, diag=TRUE)]
  
  for (block_id in blocks[upper.tri(blocks, diag=TRUE)]){
    ids=rbind(expand.grid(row_id_start[block_id]:row_id_end[block_id], 
                          col_id_start[block_id]:col_id_end[block_id]),
              expand.grid(col_id_start[block_id]:col_id_end[block_id], 
                          row_id_start[block_id]:row_id_end[block_id]))
    bigblocks[as.matrix(ids)]=block_id
  }
  
  return(bigblocks)
}


MultiBlockRescale=function(mo, pk){
  block_mean=rep(NA,length(pk))
  for (k in 1:length(pk)){
    if (k==1){
      block_mean[k]=mean(apply(mo[,1:pk[k]],2,sd))
    } else {
      if (k==length(pk)){
        block_mean[k]=mean(apply(mo[,(sum(pk[1:(k-1)])+1):ncol(mo)],2,sd))
      } else {
        block_mean[k]=mean(apply(mo[,(sum(pk[1:(k-1)])+1):sum(pk[1:(k)])],2,sd))
      }
    }
  }
  myrescaling=ifelse(block_mean/max(block_mean)<0.5, yes=1, no=0)
  if (any(myrescaling==1)){
    for (k in which(myrescaling==1)){
      if (k==1){
        mo[,1:pk[k]]=mo[,1:pk[k]]*1/(block_mean[k]/max(block_mean))
      } else {
        if (k==length(pk)){
          mo[,(sum(pk[1:(k-1)])+1):ncol(mo)]=mo[,(sum(pk[1:(k-1)])+1):ncol(mo)]*1/(block_mean[k]/max(block_mean))
        } else {
          mo[,(sum(pk[1:(k-1)])+1):sum(pk[1:k])]=mo[,(sum(pk[1:(k-1)])+1):sum(pk[1:k])]*1/(block_mean[k]/max(block_mean))
        }
      }
    }
  }
  return(mo)
}


ComputePFER=function(q, pi, N, K, method="MB"){
  if (!method%in%c("MB", "SS")){
    stop('PFER method must be "MB" or "SS".')
  }
  if (method=="MB"){
    upperbound=1/(2*pi-1)*q^2/N
  } 
  if (method=="SS"){
    ### Adapted code from stabsel package (to avoid "out of bounds" error)
    cutoff=pi
    B=K
    theta=q/N
    if (cutoff <= 3/4) {
      tmp=2 * (2 * cutoff - 1 - 1/(2*B))
    } else {
      tmp=(1 + 1/B)/(4 * (1 - cutoff + 1 / (2*B)))
    }
    upperbound=q^2/N/tmp
    if ((cutoff < 1/2 + min(theta^2, 1 / (2*B) + 3/4 * theta^2))|(cutoff>1)){
      # If out of bounds for SS formula under unimodality assumption, set to Inf
      upperbound=Inf
    }
    # upperbound=stabsel_parameters(p=N, q=q, cutoff=pi, B=K, sampling.type="SS", assumption="unimodal")$PFER
  }
  return(upperbound)
}


ComputeFDP=function(PFER, selection_proportions, pi){
  # Computing the proportion of false discoveries among discoveries
  if (length(dim(selection_proportions))==2){
    selprops=selection_proportions[upper.tri(selection_proportions)] # vector of selection proportions
  } else {
    selprops=selection_proportions
  }
  S=sum(selprops>=pi,na.rm=TRUE) # number of stability-selected edges
  if (S!=0){
    FDP=PFER/S
  } else {
    FDP=0
  }
  return(FDP)
}


GetBinomialProbabilities=function(q, N, pi, K){
  p_1=pbinom(round(K*pi), size=K, prob=q/N, log.p=TRUE) # proportion <= (1-pi)
  p_2=log(pbinom(round(K*(1-pi))-1, size=K, prob=q/N)-pbinom(round(K*pi), size=K, prob=q/N)) # 1-pi < proportion < pi
  if (is.infinite(p_2)|is.na(p_2)){
    p_2=0
    for (i in seq(round(K*(1-pi))-1, round(K*pi)+1)){
      p_2=p_2+dbinom(i, size=K, prob=q/N)
    }
    p_2=log(p_2)
  }
  p_3=pbinom(round(K*(1-pi))-1, size=K, prob=q/N, lower.tail=FALSE, log.p=TRUE) # proportion >= pi 
  
  if (abs(exp(p_1)+exp(p_2)+exp(p_3)-1)>1e-3){
    print(paste("N:", N))
    print(paste("q:", q))
    print(paste("K:", K))
    print(paste("pi:", pi))
    stop(paste0("Probabilities do not sum to 1 (Binomial distribution) \n p_1+p_2+p_3=", exp(p_1)+exp(p_2)+exp(p_3)))
  }
  return(list(p_1=p_1, p_2=p_2, p_3=p_3))
}


GetBinomialScore=function(stab_iter, q=NULL, N=NULL, pi, K){
  if (is.matrix(stab_iter)){
    stab_iter=stab_iter[upper.tri(stab_iter)]
  }
  N=length(stab_iter)
  
  if (is.null(q)){
    q=round(sum(stab_iter))
  }
  
  p_vect=GetBinomialProbabilities(q, N, 1-pi, K)
  
  S_0=sum(stab_iter<=(1-pi)) # stable-out
  S_1=sum(stab_iter>=pi) # stable-in
  U=sum((stab_iter<pi)&(stab_iter>(1-pi))) # unstable
  
  if (S_0+S_1+U!=N){
    stop(paste0("Inconsistency in number of edges \n S_0+S_1+U=", S_0+S_1+U, " instead of ", N))
  }
  
  l=S_0*p_vect$p_1+U*p_vect$p_2+S_1*p_vect$p_3
  
  if (is.infinite(l)){
    l=NA
  }
  
  return(l)
}


GetSelectionProportionsNetwork=function(data, lambda, tau=0.5, K=100, seed=1,
                                        method="graphical_lasso", implementation="glassoFast", scale=TRUE,
                                        resampling_method="subsampling", PFER_method="MB"){
  # Only method so far is graphical_lasso
  
  if ((PFER_method=="SS")&(resampling_method=="bootstrap")){
    stop("If PFER_method is set to 'SS', resampling is based on complementary pairs. It cannot be executed with the 'boostrap'. Please set resampling_method to 'subsampling' if you want to use PFER_method='SS'.")
  }
  
  set.seed(seed)
  stab_iter=matrix(0,ncol=ncol(data),nrow=ncol(data))
  
  if (!implementation%in%c("glasso", "glassoFast", "QUIC")){
    stop("For method='graphicallasso', implementation must be one of 'glasso', 'glassoFast', 'QUIC'.")
  }
  
  if (PFER_method=="MB"){
    for (i in 1:K){
      if (resampling_method=="subsampling"){
        s=sample(1:nrow(data), size=tau*nrow(data))
      } 
      if (resampling_method=="bootstrap"){
        s=sample(1:nrow(data), size=nrow(data), replace=TRUE)
      }
      data_sub=data[s,]
      
      # Making sure none of the variables has a null standard deviation
      mysd=apply(data_sub,2,sd)
      if (any(mysd==0)){
          for (k in which(mysd==0)){
              data_sub[,k]=data_sub[,k]+rnorm(n=nrow(x))
          }
      }
      
      if (scale){
        cov_sub=cor(data_sub)
      } else {
        cov_sub=cov(data_sub)
      }
      if (implementation=="glasso"){
        g_sub=glasso(s=cov_sub, rho=lambda)
        omega=g_sub$wi
      }
      if (implementation=="glassoFast"){
        g_sub=glassoFast(S=cov_sub, rho=lambda)
        omega=g_sub$wi
      }
      if (implementation=="QUIC"){ 
        g_sub=QUIC(S=cov_sub, rho=lambda, msg=0)
        omega=g_sub$X
      }
      A=ifelse(omega!=0, yes=1, no=0)
      A=A+t(A)
      A=ifelse(A!=0, yes=1, no=0)
      stab_iter=stab_iter+A
    }
    stab_iter=stab_iter/K
    diag(stab_iter)=0
  }
  
  if (PFER_method=="SS"){ # complementary pairs
    tau=0.5
    for (i in 1:ceiling(K/2)){
      # Sample 1
      s=sample(1:nrow(data), size=tau*nrow(data))
      data_sub=data[s,]
      if (scale){
        cov_sub=cor(data_sub)
      } else {
        cov_sub=cov(data_sub)
      }
      if (implementation=="glasso"){
        g_sub=glasso(s=cov_sub, rho=lambda)
        omega=g_sub$wi
      }
      if (implementation=="glassoFast"){
        g_sub=glassoFast(S=cov_sub, rho=lambda)
        omega=g_sub$wi
      }
      if (implementation=="QUIC"){ 
        g_sub=QUIC(S=cov_sub, rho=lambda, msg=0)
        omega=g_sub$X
      }
      A=ifelse(omega!=0, yes=1, no=0)
      A=A+t(A)
      A=ifelse(A!=0, yes=1, no=0)
      stab_iter=stab_iter+A
      
      # Sample 2: everything not in sample 1
      s=seq(1:nrow(data))[!seq(1:nrow(data))%in%s]
      data_sub=data[s,]
      if (scale){
        cov_sub=cor(data_sub)
      } else {
        cov_sub=cov(data_sub)
      }
      if (implementation=="glasso"){
        g_sub=glasso(s=cov_sub, rho=lambda)
        omega=g_sub$wi
      }
      if (implementation=="glassoFast"){
        g_sub=glassoFast(S=cov_sub, rho=lambda)
        omega=g_sub$wi
      }
      if (implementation=="QUIC"){ 
        g_sub=QUIC(S=cov_sub, rho=lambda, msg=0)
        omega=g_sub$X
      }
      A=ifelse(omega!=0, yes=1, no=0)
      A=A+t(A)
      A=ifelse(A!=0, yes=1, no=0)
      stab_iter=stab_iter+A
    }
    stab_iter=stab_iter/K
    diag(stab_iter)=0
  }
  
  return(stab_iter)
}


CalibrateInformationTheory=function(X, Lambda, scale=TRUE){
  # Prepare required objects
  n=nrow(X)
  if (scale){
    S=cor(X)
  } else {
    S=cov(X)
  }
  AIC=BIC=path=NULL
  
  # Loop over lambda grid
  pb=txtProgressBar(style=3)
  for (k in 1:length(Lambda)){
    setTxtProgressBar(pb, k/length(Lambda))
    myglasso=glasso(S, rho=Lambda[k])
    Theta=myglasso$wi
    loglik <- (n/2) * (log(det(Theta)) - sum(diag(Theta %*% S)))
    df <- sum(abs(Theta) > 0) - sum(abs(diag(Theta)) > 0)
    BIC=c(BIC, loglik - 0.5 * df * log(n))
    AIC=c(AIC, loglik - df)
    A=ifelse(myglasso$wi!=0, yes=1, no=0)
    A=A+t(A)
    A=ifelse(A!=0, yes=1, no=0)
    path=abind(path, A, along=3)
  }
  
  return(list(path=path, AIC=AIC, BIC=BIC))
}


StabilityCalibNetwork=function(data, pk=NULL, Lambda, lambda_other_blocks=NULL, pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, 
                               method="graphical_lasso", implementation="glassoFast", scale=TRUE,
                               resampling_method="subsampling", PFER_method="MB", PFER_thr=+Inf, FDP_thr=+Inf, 
                               verbose=TRUE, debug=FALSE){
  # Browser mode for debugging
  if (debug){
    browser()
  }
  
  # Prepare dimensions
  if (is.null(pk)){
    pk=ncol(data)
  }
  
  # Marginal correlation to get sign of the relationship
  mycor_for_sign=cor(data)
  
  # Checking consistency between arguments 
  if (!is.null(pk)){
    if (sum(pk)!=ncol(data)){
      stop("Argument pk is not consistent with the number of variables in data. Please make sure that sum(pk) is equal to ncol(data).")
    }
  }
  
  # Create matrix with block indices
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  #bigblocks_vect=factor(bigblocks_vect, levels=unique(as.vector(bigblocks))) # XXX For groups with one variable
  N_blocks=unname(table(bigblocks_vect))
  blocks=unique(as.vector(bigblocks_vect))
  names(N_blocks)=blocks
  nblocks=max(blocks)
  
  # Re-format Lambda
  if (is.vector(Lambda)){
    grid=MakeBlockLambdaGrid(Lambda=Lambda, lambda_other_blocks=lambda_other_blocks)
    Lambda=grid$Lambda
    Sequential_template=grid$Sequential_template
  } else {
    grid=MakeBlockLambdaGrid(Lambda=Lambda, lambda_other_blocks=lambda_other_blocks)
    Lambda=grid$Lambda
    Sequential_template=grid$Sequential_template
  }
  if (verbose){
    print(Lambda)
  }
  
  # Prepare the PFER and FDP thresholds
  if (length(PFER_thr)==1){
    PFER_thr_blocks=ceiling(prop.table(N_blocks)*PFER_thr)
  } else {
    if (length(PFER_thr)==nblocks){
      PFER_thr_blocks=PFER_thr
    } else {
      stop(paste0("Please provide a single number or as many values as there are blocks in PFER_thr (currently ",length(PFER_thr)," values provided)."))
    }
  }
  if (length(FDP_thr)==1){
    FDP_thr_blocks=rep(FDP_thr, nblocks)
  } else {
    if (length(FDP_thr)==nblocks){
      FDP_thr_blocks=FDP_thr
    } else {
      stop(paste0("Please provide a single number or as many values as there are blocks in FDP_thr (currently ",length(FDP_thr)," values provided)."))
    }
  }
  
  # Initialise objects to be filled
  Q=Q_s=P=matrix(NA,nrow=nrow(Lambda),ncol=nblocks)
  bigstab=array(NA,dim=c(ncol(data),ncol(data),nrow(Q)),dimnames=list(colnames(data),colnames(data),NULL))
  best_loglik=best_PFER=best_FDP=matrix(NA, nrow=nrow(Lambda), ncol=nblocks)
  if (nblocks==1){
    loglik=PFER=FDP=matrix(NA, ncol=length(pi_list), nrow=nrow(Lambda))
  } else {
    loglik=array(NA, dim=c(nrow(Lambda), length(pi_list), nblocks))
  }
  
  # Printing message 
  if (verbose){
    if (all(!is.infinite(PFER_thr_blocks))){
      print("Threshold(s) in PFER:")
      print(PFER_thr_blocks)
    }
    if (all(!is.infinite(FDP_thr_blocks))){
      print("Threshold(s) in FDP:")
      print(FDP_thr_blocks)
    }
  }
  
  # Iterate over values of lambda
  k=min_PFER=0
  if (K>1){
    if (verbose){
      pb=txtProgressBar(style=3)
    }
    while ((k<nrow(Lambda))&(min_PFER<=PFER_thr[1])&(min_PFER<=N_blocks[1])){ # 3rd entry is (very loose) stopping criterion for the FDP
      k=k+1
      if (verbose){
        setTxtProgressBar(pb, k/nrow(Lambda))
      }
      
      # Create penalisation matrix
      lambdamat=bigblocks
      for (b in 1:nblocks){
        lambdamat[bigblocks==b]=Lambda[k,b]
      }
      
      # Compute selection proportions
      stab_iter=GetSelectionProportionsNetwork(data=data, lambda=lambdamat, tau=tau, K=K, seed=seed, 
                                               method=method, implementation=implementation, scale=scale,
                                               resampling_method=resampling_method, PFER_method=PFER_method)
      bigstab[,,k]=stab_iter
      
      # Compute stability score with block-specific pi
      for (block_id in 1:nblocks){
        if (Sequential_template[k,block_id]){
          stab_iter_block=stab_iter[(bigblocks==block_id)&(upper.tri(bigblocks))] # selection proportions in the block
          q_block=round(sum(stab_iter_block)) # average number of edges selected by the original procedure in the block
          Q[k,block_id]=q_block
          N_block=length(stab_iter_block) # maximum number of edges in the block
          tmp_loglik=tmp_PFERs=tmp_FDPs=rep(NA, length(pi_list))
          
          # Compute error rates and stability score for different values of pi
          for (j in 1:length(pi_list)){
            # print(pi_list[j])
            pi=pi_list[j]
            tmp_PFERs[j]=ComputePFER(q=q_block, pi=pi, N=N_block, K=K, method=PFER_method)
            tmp_FDPs[j]=ComputeFDP(PFER=tmp_PFERs[j], pi=pi, selection_proportions=stab_iter_block)
            if ((tmp_PFERs[j]<=PFER_thr_blocks[block_id])&(tmp_FDPs[j]<FDP_thr_blocks[block_id])){
              tmp_loglik[j]=GetBinomialScore(stab_iter=stab_iter_block, q=q_block, N=N_block, pi=pi, K=K)
              # print(tmp_loglik[j])
            }
          }
          
          # Store stability score in a matrix if only one block
          if (nblocks==1){
            loglik[k,]=tmp_loglik
            PFER[k,]=tmp_PFERs
            FDP[k,]=tmp_FDPs
          } else {
            loglik[k,,block_id]=tmp_loglik
          }
          
          # Keep best stability score and other parameters at the max
          if (any(!is.na(tmp_loglik))){
            tmp_loglik[is.na(tmp_loglik)]=0
            myid=which.min(tmp_loglik)
            tmp_loglik[which(tmp_loglik==0)]=NA
            best_loglik[k,block_id]=tmp_loglik[myid]
            P[k,block_id]=pi_list[myid]
            Q_s[k,block_id]=sum(stab_iter_block>=pi_list[myid])
            best_PFER[k,block_id]=tmp_PFERs[myid]
            best_FDP[k,block_id]=tmp_FDPs[myid]
            if (nblocks==1){ # Stopping criteria are used only for single-block calibration
              min_PFER=tmp_PFERs[myid]
            }
          }
        }
      }
    } 
    if (verbose){
      setTxtProgressBar(pb, 1)
      cat("\n")
    }
    best_loglik_blocks=best_loglik
    best_loglik=matrix(apply(best_loglik,1,sum), ncol=1)
  } else {
    PFER=matrix(NA,nrow(Lambda),length(pi_list))
    for (k in 1:nrow(Lambda)){
      # Create penalisation matrix
      lambdamat=bigblocks
      for (b in 1:nblocks){
        lambdamat[bigblocks==b]=Lambda[k,b]
      }
      
      # Compute selection proportions
      stab_iter=GetSelectionProportionsNetwork(data=data, lambda=lambdamat, tau=tau, K=K, seed=seed, 
                                               method=method, implementation=implementation, scale=scale,
                                               resampling_method=resampling_method, PFER_method=PFER_method)
      
      # Get number of selected variables per block
      for (block_id in 1:nblocks){
        stab_iter_block=stab_iter[(bigblocks==block_id)&(upper.tri(bigblocks))] # selection proportions in the block
        q_block=round(sum(stab_iter_block)) # average number of edges selected by the original procedure in the block
        Q[k,block_id]=q_block
      }
    }
  }
  
  # Prepare outputs
  if (K>1){
    if (nblocks==1){
      return(list(S=-best_loglik, Lambda=Lambda, 
                  Q=Q, Q_s=Q_s, P=P,
                  PFER=best_PFER, FDP=best_FDP,
                  S_2d=-loglik, PFER_2d=PFER, FDP_2d=FDP, 
                  selprop=bigstab,
                  sign=sign(mycor_for_sign),
                  methods=list(method=method, implementation=implementation, resampling_method=resampling_method, PFER_method=PFER_method),
                  params=list(K=K, pi_list=pi_list, tau=tau, pk=pk, PFER_thr=PFER_thr, FDP_thr=FDP_thr, seed=seed, data=data)))
    } else {
      return(list(S_blocks=-best_loglik_blocks, Lambda=Lambda, 
                  Q=Q, Q_s=Q_s, P=P,
                  PFER=best_PFER, FDP=best_FDP,
                  S_2d=-loglik,
                  selprop=bigstab,
                  sign=sign(mycor_for_sign),
                  methods=list(method=method, implementation=implementation, resampling_method=resampling_method, PFER_method=PFER_method),
                  params=list(K=K, pi_list=pi_list, tau=tau, pk=pk, PFER_thr=PFER_thr, FDP_thr=FDP_thr, seed=seed, data=data)))
    }
  } else {
    return(list(Q=Q))
  }
}


StabilityCalibRegression=function(xdata, ydata, pk=NULL, Lambda, pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, 
                                  family="gaussian", implementation="glmnet",
                                  resampling_method="subsampling", resampling_function=NULL, PFER_method="MB", PFER_thr=+Inf, FDP_thr=+Inf, 
                                  verbose=TRUE, debug=FALSE, ...){
  # Browser mode for debugging
  if (debug){
    browser()
  }
  
  # Prepare dimensions
  if (is.null(pk)){
    pk=ncol(xdata)
  }
  
  # Checking consistency between arguments 
  if (!is.null(pk)){
    if (sum(pk)!=ncol(xdata)){
      stop("Argument pk is not consistent with the number of variables in X Please make sure that sum(pk) is equal to ncol(X).")
    }
  }
  
  # Refine K if using complementary pairs (SS)
  if (PFER_method=="SS"){
    K=ceiling(K/2)*2
  }
  
  # Name columns of xdata
  if (is.null(colnames(xdata))){
    colnames(xdata)=paste0("var", 1:ncol(xdata))
  }
  
  # Create matrix with block indices
  bigblocks_vect=diag(GetBlockMatrix(pk))
  N_blocks=unname(table(bigblocks_vect))
  blocks=unique(bigblocks_vect)
  names(N_blocks)=blocks
  nblocks=max(blocks)
  
  # Re-format Lambda
  if (is.vector(Lambda)){
    Lambda=matrix(rep(Lambda, nblocks), ncol=nblocks)
  }
  rownames(Lambda)=paste0("s",seq(0,nrow(Lambda)-1))
  
  # Re-format ydata
  if (is.vector(ydata)|is.factor(ydata)){
    ydata=matrix(ydata, ncol=1)
  }
  
  # Check consistency between X and Y
  if (nrow(xdata)!=nrow(ydata)){
    stop("Different numbers of observations in xdata and ydata.")
  }
  
  # Prepare the PFER and FDP thresholds
  if (length(PFER_thr)==1){
    PFER_thr_blocks=ceiling(prop.table(N_blocks)*PFER_thr)
  } else {
    if (length(PFER_thr)==nblocks){
      PFER_thr_blocks=PFER_thr
    } else {
      stop(paste0("Please provide a single number or as many values as there are blocks in PFER_thr (currently ",length(PFER_thr)," values provided)."))
    }
  }
  if (length(FDP_thr)==1){
    FDP_thr_blocks=rep(FDP_thr, nblocks)
  } else {
    if (length(FDP_thr)==nblocks){
      FDP_thr_blocks=FDP_thr
    } else {
      stop(paste0("Please provide a single number or as many values as there are blocks in FDP_thr (currently ",length(FDP_thr)," values provided)."))
    }
  }
  
  # Initialise objects to be filled
  Q=P=Q_s=matrix(NA, nrow=nrow(Lambda), ncol=1) # average number of included variable per lambda
  loglik=PFER=FDP=matrix(NA, nrow=nrow(Lambda), ncol=length(pi_list))
  N=N_block=ncol(xdata)
  Beta=array(0, dim=c(nrow(Lambda), ncol(xdata), K))
  rownames(Beta)=rownames(Lambda)
  colnames(Beta)=colnames(xdata)
  if (family=="mgaussian"){
    Beta_mgaussian=array(0, dim=c(nrow(Lambda), ncol(xdata), K, ncol(ydata)))
    rownames(Beta_mgaussian)=rownames(Lambda)
    colnames(Beta_mgaussian)=colnames(xdata)
  }
  best_loglik=best_PFER=best_FDP=matrix(NA, nrow=nrow(Lambda), ncol=nblocks)
  
  # Printing message 
  if (verbose){
    if (all(!is.infinite(PFER_thr_blocks))){
      print("Threshold(s) in PFER:")
      print(PFER_thr_blocks)
    }
    if (all(!is.infinite(FDP_thr_blocks))){
      print("Threshold(s) in FDP:")
      print(FDP_thr_blocks)
    }
  }
  
  # Computation of the selection proportions over Lambda
  pb=txtProgressBar(style=3)
  if (PFER_method=="MB"){
    for (k in 1:K){
      # print(k)
      setTxtProgressBar(pb, k/K)
      set.seed(k)
      s=GetSubsample(ydata=ydata, family=family, tau=tau, resampling_method=resampling_method, resampling_function=resampling_function)
      Xsub = xdata[s,]
      Ysub = ydata[s,]
      # model_sub=glmnet(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, ...)
      # if (family=="mgaussian"){
      #   beta_sub=coef(model_sub)[[1]]
      # } else {
      #   beta_sub=coef(model_sub)
      # }
      # beta_sub=t(as.matrix(beta_sub))
      # beta_sub=beta_sub[,colnames(xdata)] # removing the intercept if included
      # if (family=="mgaussian"){
      #   beta_sub=ifelse(beta_sub!=0,yes=1,no=0)
      # }
      # Beta[rownames(beta_sub),colnames(beta_sub),k]=beta_sub
      # if (family=="mgaussian"){
      #   for (y_id in 1:ncol(ydata)){
      #     beta_sub=coef(model_sub)[[y_id]]
      #     beta_sub=t(as.matrix(beta_sub))
      #     beta_sub=beta_sub[,colnames(xdata)] # removing the intercept if included
      #     Beta_mgaussian[rownames(beta_sub),colnames(beta_sub),k,y_id]=beta_sub
      #   }
      # }
      
      mybeta=SelectionFunction(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)
      if (family=="mgaussian"){
        beta_sub=mybeta[,,1]
        beta_sub=ifelse(beta_sub!=0,yes=1,no=0)
      } else {
        beta_sub=mybeta
      }
      Beta[rownames(beta_sub),colnames(beta_sub),k]=beta_sub
      if (family=="mgaussian"){
        for (y_id in 1:ncol(ydata)){
          beta_sub=mybeta[,,y_id]
          Beta_mgaussian[rownames(beta_sub),colnames(beta_sub),k,y_id]=beta_sub
        }
      }
    }
  }
  if (PFER_method=="SS"){
    for (k in 1:ceiling(K/2)){
      setTxtProgressBar(pb, 2*k/K)
      set.seed(k)
      s=GetSubsample(ydata=ydata, family=family, tau=tau, resampling_method=resampling_method, resampling_function=resampling_function)
      
      # Fist subset
      Xsub = xdata[s,]
      Ysub = ydata[s,]
      mybeta=SelectionFunction(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)
      if (family=="mgaussian"){
        beta_sub=mybeta[,,1]
        beta_sub=ifelse(beta_sub!=0,yes=1,no=0)
      } else {
        beta_sub=mybeta
      }
      Beta[rownames(beta_sub),colnames(beta_sub),k]=beta_sub
      if (family=="mgaussian"){
        for (y_id in 1:ncol(ydata)){
          beta_sub=mybeta[,,y_id]
          Beta_mgaussian[rownames(beta_sub),colnames(beta_sub),k,y_id]=beta_sub
        }
      }
      
      # Complementary subset
      Xsub = xdata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
      Ysub = ydata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
      mybeta=SelectionFunction(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)
      if (family=="mgaussian"){
        beta_sub=mybeta[,,1]
        beta_sub=ifelse(beta_sub!=0,yes=1,no=0)
      } else {
        beta_sub=mybeta
      }
      Beta[rownames(beta_sub),colnames(beta_sub),ceiling(K/2)+k]=beta_sub
      if (family=="mgaussian"){
        for (y_id in 1:ncol(ydata)){
          beta_sub=mybeta[,,y_id]
          Beta_mgaussian[rownames(beta_sub),colnames(beta_sub),ceiling(K/2)+k,y_id]=beta_sub
        }
      }
    }
  }
  
  cat("\n")
  if (K>1){
    bigstab=apply(Beta,c(1,2),FUN=function(x){sum(x!=0)})/K # selection proportions
  }
  
  # Computation of the stability score over Lambda and pi_list
  if (K>1){
    for (k in 1:nrow(Lambda)){
      q_block=mean(apply(Beta[k,,],2,FUN=function(x){sum(x!=0)}))
      Q[k,1]=round(q_block)
      for (j in 1:length(pi_list)){
        pi=pi_list[j]
        PFER[k,j]=ComputePFER(q=q_block, pi=pi, N=N_block, K=K, method=PFER_method)
        FDP[k,j]=ComputeFDP(PFER=PFER[k,j], pi=pi, selection_proportions=bigstab[k,])
        if ((PFER[k,j]<=PFER_thr)&(FDP[k,j]<=FDP_thr)){
          loglik[k,j]=GetBinomialScore(stab_iter=bigstab[k,], pi=pi, K=K)
        }
      }
    }
    rownames(loglik)=rownames(PFER)=rownames(FDP)=rownames(bigstab)
    colnames(loglik)=colnames(PFER)=colnames(FDP)=pi_list
    
    # Add constraint
    if ((!is.infinite(PFER_thr))|(!is.infinite(FDP_thr))){
      if (!is.infinite(PFER_thr)){
        loglik=ifelse(PFER>PFER_thr, yes=NA, no=loglik)
      } else {
        loglik=ifelse(FDP>FDP_thr, yes=NA, no=loglik)
      }
    }
    
    # Computation of the best score by lambda
    block_id=1
    for (k in 1:nrow(Lambda)){
      tmp_loglik=loglik[k,]
      if (any(!is.na(tmp_loglik))){
        tmp_loglik[is.na(tmp_loglik)]=0
        myid=which.min(tmp_loglik)
        tmp_loglik[which(tmp_loglik==0)]=NA
        best_loglik[k,block_id]=tmp_loglik[myid]
        P[k,block_id]=pi_list[myid]
        Q_s[k,block_id]=sum(bigstab[k,]>=pi_list[myid])
        best_PFER[k,block_id]=PFER[k,myid]
        best_FDP[k,block_id]=FDP[k,myid]
        # if (nblocks==1){ # Stopping criteria are not used only for single-block calibration
        #   min_PFER=PFER[k,myid]
        # }
      }
    }
  } else {
    for (k in 1:nrow(Lambda)){
      q_block=sum(myscreen$Beta[k,,1]!=0)
      Q[k,1]=round(q_block)
      for (j in 1:length(pi_list)){
        pi=pi_list[j]
        PFER[k,j]=ComputePFER(q=q_block,pi=pi,N=N_block)
      }
    }
  }
  
  if (family=="mgaussian"){
    Beta=Beta_mgaussian
  }
  
  # Prepare outputs
  if (nblocks==1){ # Only possible with 1 block so far
    if (K>1){
      return(list(S=-best_loglik, Lambda=Lambda, 
                  Q=Q, Q_s=Q_s, P=P,
                  PFER=best_PFER, FDP=best_FDP,
                  S_2d=-loglik, PFER_2d=PFER, FDP_2d=FDP, 
                  selprop=bigstab, Beta=Beta,
                  methods=list(family=family, resampling_method=resampling_method, PFER_method=PFER_method),
                  params=list(K=K, pi_list=pi_list, tau=tau, pk=pk, PFER_thr=PFER_thr, FDP_thr=FDP_thr, seed=seed, xdata=xdata, ydata=ydata)))
    } else {
      return(list(Q=Q, Beta=Beta, PFER_2d=PFER))
    }
  }
}


LambdaGridNetwork=function(data, pk=NULL, lambda_other_blocks=NULL, pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, 
                           method="graphical_lasso", implementation="glassoFast", scale=TRUE,
                           resampling_method="subsampling", PFER_method="MB", PFER_thr=+Inf, FDP_thr=+Inf, 
                           lambda_path_factor=0.0001, max_density=0.3, 
                           lambda_path_refined_cardinal=50, verbose=TRUE){
  # Ideally for multi-block calibration, set lambda_other_blocks to a null vector of length the number of blocks
  # However, very long to compute for large number of nodes 
  # Instead: approximation where the lambda_dense for other blocks is one that gives density 0.3 (or above given PFER_thr or FDP_thr)
  # Working well when the network is not expected to be too dense 
  
  # Prepare dimensions
  if (is.null(pk)){
    pk=ncol(data)
  }
  p=sum(pk)
  N=p*(p-1)/2
  
  # Prepare lambda_dense
  ldense=lambda_other_blocks
  
  # Create matrix with block indices
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  N_blocks=unname(table(bigblocks_vect))
  blocks=unique(as.vector(bigblocks_vect))
  names(N_blocks)=blocks
  nblocks=max(blocks)
  
  # Prepare the PFER and FDP thresholds
  if (length(PFER_thr)==1){
    PFER_thr_blocks=ceiling(prop.table(N_blocks)*PFER_thr)
  } else {
    if (length(PFER_thr)==nblocks){
      PFER_thr_blocks=PFER_thr
    } else {
      stop(paste0("Please provide a single number or as many values as there are blocks in PFER_thr (currently ",length(PFER_thr)," values provided)."))
    }
  }
  if (length(FDP_thr)==1){
    FDP_thr_blocks=rep(FDP_thr, nblocks)
  } else {
    if (length(FDP_thr)==nblocks){
      FDP_thr_blocks=FDP_thr
    } else {
      stop(paste0("Please provide a single number or as many values as there are blocks in FDP_thr (currently ",length(FDP_thr)," values provided)."))
    }
  }
  
  # Printing message 
  if (verbose){
    if (all(!is.infinite(PFER_thr_blocks))){
      print("Threshold(s) in PFER:")
      print(PFER_thr_blocks)
    }
    if (all(!is.infinite(FDP_thr_blocks))){
      print("Threshold(s) in FDP:")
      print(FDP_thr_blocks)
    }
  }
  
  # Get upperbound of Lambda
  if (scale){
    mycov=cor(data)
  } else {
    mycov=cov(data)
  }
  lmax=getMaxCov(mycov)
  lmin=lmax
  lmin=rep(lmin, nblocks)
  
  # Identifying the constraint
  if (all(is.infinite(PFER_thr))){
    if (all(is.infinite(FDP_thr))){
      type_opt_problem="unconstrained"
    } else {
      type_opt_problem="constrained_FDP"
    }
  } else {
    type_opt_problem="constrained_PFER"
  }
  
  if (type_opt_problem=="unconstrained"){
    max_q=rep(0,nblocks)
    redo=TRUE
    done=rep(0, nblocks)
    while (redo){
      lmin=lmin*lambda_path_factor
      Lambda=NULL
      for (b in 1:nblocks){
        Lambda=cbind(Lambda, GetLambdaPath(lmax, lmin[b], cardinal=lambda_path_refined_cardinal))
      }
      
      # Initialisation of the smallest lambda
      lmin=Lambda[2,]
      l=1
      while (l<nrow(Lambda)){
        if (is.null(lambda_other_blocks)){
          ldense=lmin
        }
        tmpLambda=Lambda[c(1,l),]
        myscreen=StabilityCalibNetwork(data=data, pk=pk, Lambda=tmpLambda, lambda_other_blocks=ldense, pi_list=pi_list, K=1, tau=tau, seed=seed, 
                                       method=method, implementation=implementation, resampling_method=resampling_method, scale=scale,
                                       PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                                       verbose=FALSE, debug=FALSE) # Only 1 iteration to get the Q
        
        if (l<nrow(Lambda)){
          # Updating the smallest lambda if the density of the block is still below max_density
          # lmin=ifelse(myscreen$Q[2,]<=max_density*N_blocks, yes=Lambda[l+1,], no=lmin)
          for (b in 1:nblocks){
            lmin[b]=ifelse((myscreen$Q[seq(2,nrow(myscreen$Q),by=2)[b],b]<=(max_density*N_blocks)[b])&(done[b]==0), 
                           yes=Lambda[l+1,b], no=lmin[b])
            done[b]=ifelse(myscreen$Q[seq(2,nrow(myscreen$Q),by=2)[b],b]>(max_density*N_blocks)[b], yes=1, no=0)
          }
        }
        
        # Increment if max_density is not yet reached
        Q_block_iteration=NULL
        for (b in 1:nblocks){
          Q_block_iteration=c(Q_block_iteration, myscreen$Q[seq(2,nrow(myscreen$Q),by=2)[b],b])
        }
        # if (all(is.na(Q_block_iteration))){
        #   Q_block_iteration=0
        # }
        # print(Q_block_iteration)
        if (any(Q_block_iteration<=max_density*N_blocks)){
          l=l+1
        } else {
          l=nrow(Lambda) # stopping current while loop 
          redo=FALSE # stopping overarching while loop
        }
      }
    }
  }
  
  if (type_opt_problem=="constrained_PFER"){
    max_q=rep(0,nblocks)
    redo=TRUE
    print(PFER_thr)
    done=rep(0, nblocks)
    while (redo){
      # print(paste0("lmin: ", lmin))
      lmin=lmin*lambda_path_factor
      Lambda=NULL
      for (b in 1:nblocks){
        Lambda=cbind(Lambda, GetLambdaPath(lmax, lmin[b], cardinal=lambda_path_refined_cardinal))
      }
      # print("Lambda:")
      # print(Lambda)
      
      # Initialisation of the smallest lambda
      lmin=Lambda[2,]
      l=1
      while (l<nrow(Lambda)){
        if (is.null(lambda_other_blocks)){
          ldense=lmin
        }
        tmpLambda=Lambda[c(1,l),,drop=FALSE]
        tmpLambda[2,]=ifelse(tmpLambda[2,]<lmin, yes=lmin, no=tmpLambda[2,])
        myscreen=StabilityCalibNetwork(data=data, pk=pk, Lambda=tmpLambda, lambda_other_blocks=ldense, pi_list=pi_list, K=1, tau=tau, seed=seed, 
                                       method=method, implementation=implementation, resampling_method=resampling_method, scale=scale,
                                       PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                                       verbose=FALSE, debug=FALSE) # Only 1 iteration to get the Q
        
        # Compute PFER
        PFER_l=rep(NA, nblocks)
        for (b in 1:nblocks){
          mytmplist=NULL
          for (j in 1:length(pi_list)){
            pi=pi_list[j]
            mytmplist=c(mytmplist, ComputePFER(q=myscreen$Q[2,b], pi=pi, N=N_blocks[b], K=K, method=PFER_method))
          }
          PFER_l[b]=min(mytmplist)
        }
        # print(l)
        # print(print(nrow(Lambda)))
        # print("---")
        # print(paste0("Lambda: ", tmpLambda[2,]))
        # print(paste0("Q: ", myscreen$Q[2,]))
        # print(paste0("PFER: ", PFER_l))
        
        if ((l<nrow(Lambda))&(is.null(lambda_other_blocks))){
          # Updating the smallest lambda if the PFER of the block is still below the threshold (with some margin)
          lmin=ifelse((PFER_l<=(PFER_thr_blocks*1.2+1))&(done==0), yes=Lambda[l+1,], no=lmin)
          done=ifelse(PFER_l>(PFER_thr_blocks*1.2+1), yes=1, no=0)
        }
        # if (l<nrow(Lambda)){
        #   # Updating the smallest lambda if the density of the block is still below max_density
        #   # lmin=ifelse(myscreen$Q[2,]<=max_density*N_blocks, yes=Lambda[l+1,], no=lmin)
        #   for (b in 1:nblocks){
        #     lmin[b]=ifelse(myscreen$Q[seq(2,nrow(myscreen$Q),by=2)[b],b]<=(max_density*N_blocks)[b], yes=Lambda[l+1,b], no=lmin[b])
        #   }
        # }
        # print(paste0("lmin: ", lmin))
        # print(done)
        
        # Increment if PFER or max_density are not yet reached
        Q_block_iteration=NULL
        for (b in 1:nblocks){
          Q_block_iteration=c(Q_block_iteration, myscreen$Q[seq(2,nrow(myscreen$Q),by=2)[b],b])
        }
        # print(paste0("Q_block: ", Q_block_iteration))
        # print(max_density*N_blocks)
        # cat("\n")
        
        if (any(PFER_l<=(PFER_thr_blocks*1.2+1))&(any(Q_block_iteration<=max_density*N_blocks))){
          l=l+1
        } else {
          l=nrow(Lambda) # stopping current while loop 
          redo=FALSE # stopping overarching while loop
        }
      }
    }
  }
  
  if (type_opt_problem=="constrained_FDP"){
    max_q=rep(0,nblocks)
    redo=TRUE
    done=rep(0, nblocks)
    while (redo){
      lmin=lmin*lambda_path_factor
      Lambda=NULL
      for (b in 1:nblocks){
        Lambda=cbind(Lambda, GetLambdaPath(lmax, lmin[b], cardinal=lambda_path_refined_cardinal))
      }
      
      # Initialisation of the smallest lambda
      lmin=Lambda[2,]
      l=1
      while (l<nrow(Lambda)){
        if (is.null(lambda_other_blocks)){
          ldense=lmin
        }
        tmpLambda=Lambda[c(1,l),,drop=FALSE]
        tmpLambda[2,]=ifelse(tmpLambda[2,]<lmin, yes=lmin, no=tmpLambda[2,])
        myscreen=StabilityCalibNetwork(data=data, pk=pk, Lambda=tmpLambda, lambda_other_blocks=ldense, pi_list=pi_list, K=1, tau=tau, seed=seed, 
                                       method=method, implementation=implementation, resampling_method=resampling_method, scale=scale,
                                       PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                                       verbose=FALSE, debug=FALSE) # Only 1 iteration to get the Q
        
        # Compute PFER
        PFER_l=rep(NA, nblocks)
        for (b in 1:nblocks){
          mytmplist=NULL
          for (j in 1:length(pi_list)){
            pi=pi_list[j]
            mytmplist=c(mytmplist, ComputePFER(q=myscreen$Q[2,b], pi=pi, N=N_blocks[b], K=K, method=PFER_method))
          }
          PFER_l[b]=min(mytmplist)
        }
        # print(tmpLambda[2,])
        # print(myscreen$Q[2,])
        # print(PFER_l)
        
        # if ((l<nrow(Lambda))&(is.null(lambda_other_blocks))){
        #   # Updating the smallest lambda if the PFER of the block is still below the threshold (with some margin)
        #   lmin=ifelse(PFER_l<=N, yes=Lambda[l+1,], no=lmin)
        # }
        if ((l<nrow(Lambda))&(is.null(lambda_other_blocks))){
          # Updating the smallest lambda if the PFER of the block is still below the threshold (with some margin)
          lmin=ifelse((PFER_l<=(N_blocks*1.2+1))&(done==0), yes=Lambda[l+1,], no=lmin)
          done=ifelse(PFER_l>(N_blocks*1.2+1), yes=1, no=0)
        }
        
        # Increment if max_density is not yet reached
        if (any(PFER_l<=N)){
          l=l+1
        } else {
          l=nrow(Lambda) # stopping current while loop 
          redo=FALSE # stopping overarching while loop
        }
      }
    }
  }
  
  # Prepare final lambda path for each block
  Lambda=NULL
  for (b in 1:nblocks){
    Lambda=cbind(Lambda, GetLambdaPath(lmax, lmin[b], cardinal=lambda_path_refined_cardinal))
  }
  
  if (ncol(Lambda)==1){
    Lambda=Lambda[,1]
  } 
  
  return(Lambda)
}


GetDensityQ=function(out){ # Returns block specific density with the original procedure (for networks)
  pk=out$params$pk
  M=GetBlockMatrix(pk)
  density=matrix(NA,ncol=3,nrow=dim(out$selprop)[3])
  for (k in 1:dim(out$selprop)[3]){
    for (b in 1:3){
      density[k,b]=sum(out$selprop[,,k][which((M==b)&upper.tri(M))])/sum((M==b)&upper.tri(M))
    }
  }
  return(density)
}


LambdaGridRegression=function(xdata, ydata, tau=0.5, seed=1, 
                              family="gaussian", implementation="glmnet", 
                              resampling_method="subsampling", resampling_function=NULL, PFER_method="MB", PFER_thr=+Inf, FDP_thr=+Inf, 
                              lambda_path_factor=0.0001, max_density=0.3, 
                              lambda_path_refined_cardinal=100, verbose=TRUE, ...){
  # Prepare dimensions
  pk=ncol(xdata)
  N=pk
  
  # Prepare data format
  if (is.vector(ydata)|is.factor(ydata)){
    ydata=matrix(ydata, ncol=1)
  }
  
  # Check consistency between X and Y
  if (nrow(xdata)!=nrow(ydata)){
    stop("Different numbers of observations in xdata and ydata.")
  }
  
  # Name columns of xdata
  if (is.null(colnames(xdata))){
    colnames(xdata)=paste0("var", 1:ncol(xdata))
  }
  
  # Get upperbound of Lambda
  set.seed(1)
  s=GetSubsample(ydata=ydata, family=family, tau=tau, resampling_method=resampling_method, resampling_function=resampling_function)
  set.seed(1)
  if (implementation=="glmnet"){
    mycv=cv.glmnet(x=xdata[s,], y=ydata[s,], family=family, ...)
  }
  if (implementation=="gglasso"){
    if (family=="binomial"){
      ytmp=ydata[s,]
      ytmp[ytmp==min(ytmp)]=-1
      ytmp[ytmp==max(ytmp)]=1
      mycv=cv.gglasso(x=xdata[s,], y=ytmp, loss="logit", ...)
    } 
    if (family=="gaussian"){
      mycv=cv.gglasso(x=xdata[s,], y=ydata[s,], loss="ls", ...)
    }
    if (!family%in%c("binomial", "gaussian")){
      stop("Please set family as 'binomial' or 'gaussian' for use with gglasso.")
    }
    group=mycv$gglasso.fit$group
    mycv$nzero=apply(mycv$gglasso.fit$beta,2,FUN=function(M){
      mysum=0
      for (tmpgroup in unique(group)){
        if (any(M[group==tmpgroup]!=0)){
          mysum=mysum+1
        }
      }
      return(mysum)
    })
    
  }
  if (verbose){
    plot(mycv, las=1)
    if (implementation=="glmnet"){
      print(paste("Minimum number of selected variables:", min(mycv$nzero)))
      print(paste("Maximum number of selected variables:", max(mycv$nzero)))
    }
    if (implementation=="gglasso"){
      print(paste("Minimum number of selected groups:", min(mycv$nzero)))
      print(paste("Maximum number of selected groups:", max(mycv$nzero)))
    }
  }
  Lambda=cbind(GetLambdaPath(lmax=max(mycv$lambda), lmin=min(mycv$lambda), cardinal=lambda_path_refined_cardinal))
  
  return(Lambda)
}


RefinedGridCalibNetwork=function(data, pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, 
                                 method="graphical_lasso", implementation="glassoFast", scale=TRUE,
                                 resampling_method="subsampling", PFER_method="MB", PFER_thr=+Inf, FDP_thr=+Inf, 
                                 lambda_path_factor=0.0001, max_density=0.3, 
                                 lambda_path_refined_cardinal=50, stopping_criterion="Q", 
                                 verbose=TRUE, debug=FALSE){
  # Get data dimensions
  pk=ncol(data)
  N=sum(pk)*(sum(pk)-1)/2
  
  # Checking inputs
  if (!stopping_criterion%in%c("Q", "Q_s")){
    stop("Argument stopping_criterion must be set to Q (number of edges selected by the original procedure) or Q_s (number of stable edges).")
  }
  if (length(PFER_thr)>1){
    stop("Please provide a single number for PFER_thr.")
  }
  if (length(FDP_thr)>1){
    stop("Please provide a single number for FDP_thr.")
  }
  
  # Create matrix with block indices
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  N_blocks=unname(table(bigblocks_vect))
  blocks=unique(as.vector(bigblocks_vect))
  names(N_blocks)=blocks
  nblocks=max(blocks)
  
  # Prepare the PFER and FDP thresholds
  if (length(PFER_thr)==1){
    PFER_thr_blocks=ceiling(prop.table(N_blocks)*PFER_thr)
  } else {
    if (length(PFER_thr)==nblocks){
      PFER_thr_blocks=PFER_thr
    } else {
      stop(paste0("Please provide a single number or as many values as there are blocks in PFER_thr (currently ",length(PFER_thr)," values provided)."))
    }
  }
  if (length(FDP_thr)==1){
    FDP_thr_blocks=rep(FDP_thr, nblocks)
  } else {
    if (length(FDP_thr)==nblocks){
      FDP_thr_blocks=FDP_thr
    } else {
      stop(paste0("Please provide a single number or as many values as there are blocks in FDP_thr (currently ",length(FDP_thr)," values provided)."))
    }
  }
  
  # Define a broad grid of lambda values
  if (verbose){
    print("Defining a broad grid of lambda values...")
  }
  Lambda=LambdaGridNetwork(data=data, tau=tau, seed=seed, 
                           method=method, implementation=implementation, scale=scale,
                           resampling_method=resampling_method, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                           lambda_path_factor=lambda_path_factor, max_density=max_density, 
                           lambda_path_refined_cardinal=lambda_path_refined_cardinal, verbose=FALSE)
  if (verbose){
    print("Broad grid of lambda values:")
    print(Lambda)
  }
  
  # Iterate over the lambda paths until as fine as possible
  delta=Inf
  k=0
  while (delta>1){
    k=k+1
    if (verbose){
      print(paste0("Iteration ", k, "..."))
    }
    out=StabilityCalibNetwork(data=data, pk=pk, Lambda=Lambda, lambda_other_blocks=NULL, pi_list=pi_list, K=K, tau=tau, seed=seed, 
                              method=method, implementation=implementation, resampling_method=resampling_method, scale=scale,
                              PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                              verbose=FALSE, debug=FALSE)
    if (verbose){
      print("Corresponding Q:")
      print(as.vector(out$Q))
      print("Corresponding number of stable edges:")
      print(as.vector(out$Q_s))
    }
    
    # Extract best lambda of current iteration
    myS=out$S[,1] # Best stability score 
    myS[is.na(myS)]=0
    id=which.max(myS) # ID of the best stability score around which it will be continued
    id_bottom=ifelse((id-1)>=1, yes=id-1, no=1)
    id_top=ifelse((id+1)<=length(myS), yes=id+1, no=length(myS))
    if (stopping_criterion=="Q"){
      delta=max(c(out$Q[id,1]-out$Q[id_bottom,1], out$Q[id_top,1]-out$Q[id,1]))
    } 
    if (stopping_criterion=="Q_s"){
      delta=max(c(out$Q_s[id,1]-out$Q_s[id_bottom,1], out$Q_s[id_top,1]-out$Q_s[id,1]))
    }
    
    # Prepare refined grid of lambda values for next iteration
    lmax=Lambda[id_bottom]
    lmin=Lambda[id_top]
    Lambda=GetLambdaPath(lmax,lmin,cardinal=lambda_path_refined_cardinal) 
    
    # Storing the outputs in a large object
    if (k==1){
      out_full=out
    } else {
      out_full=MergeOutputsNetwork(out_full, out)
    }
  }
  
  # Return the outputs
  # return(list(output_full=out_full, output_final=out))
  return(out_full)
}


CalibrateNetwork=function(data, pk=NULL, refine_calib_grid=TRUE,
                          Lambda=NULL, lambda_other_blocks=NULL, pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, 
                          method="graphical_lasso", implementation="glassoFast", scale=TRUE,
                          resampling_method="subsampling", PFER_method="MB", PFER_thr=+Inf, FDP_thr=+Inf, 
                          lambda_path_factor=0.0001, max_density=0.3, 
                          lambda_path_refined_cardinal=50, 
                          stopping_criterion="Q", using_empty_blocks=TRUE,
                          verbose=TRUE, debug=FALSE){
  # Prepare objects
  if (is.null(pk)){
    pk=ncol(data)
  }
  if (length(pk)>1){
    calibration="multi-block"
  } else {
    calibration="single-block"
  }
  if (verbose){
    print(paste("Starting", calibration, "calibration..."))
  }
  
  # Create matrix with block indices
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  blocks=unique(as.vector(bigblocks_vect))
  nblocks=max(blocks)
  
  # Define use of refined grid
  if ((calibration=="multi-block")|(!is.null(Lambda))){
    refine_calib_grid=FALSE
  }
  if (verbose){
    if (refine_calib_grid){
      print("Using iteratively refined grids of lambda")
    } else {
      print("Not using iteratively refined grids of lambda")
    }
  }
  
  # Checking inputs
  if (!is.null(Lambda)){
    if (is.matrix(Lambda)){
      if ((ncol(Lambda)!=nblocks)&(ncol(Lambda)!=1)){
        stop(paste0("Please provide a matrix Lambda with as many columns as blocks (N=",nblocks,")."))
      }
    }
  }
  
  # Launch calibration
  if (refine_calib_grid){
    # Single-block calibration with refined grids
    out=RefinedGridCalibNetwork(data=data, pi_list=pi_list, K=K, tau=tau, seed=seed, 
                                method=method, implementation=implementation, scale=scale,
                                resampling_method=resampling_method, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                                lambda_path_factor=lambda_path_factor, max_density=max_density, 
                                lambda_path_refined_cardinal=lambda_path_refined_cardinal, 
                                stopping_criterion=stopping_criterion, 
                                verbose=verbose, debug=FALSE)
  } else {
    if (is.null(Lambda)){
      # Define a broad grid of lambda values
      if (verbose){
        print("Defining a broad grid of lambda values...")
      }
      Lambda=LambdaGridNetwork(data=data, tau=tau, seed=seed, 
                               method=method, implementation=implementation, scale=scale,
                               resampling_method=resampling_method, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                               lambda_path_factor=lambda_path_factor, max_density=max_density, 
                               lambda_path_refined_cardinal=lambda_path_refined_cardinal, verbose=FALSE)
      if (verbose){
        print("Broad grid of lambda values:")
        print(Lambda)
      }
      
      # Re-formatting Lambda 
      Lambda=matrix(rep(Lambda, nblocks), ncol=nblocks)
      Lambda=expand.grid(data.frame(Lambda))
    } 
    
    # Single or multi-block calibration on grid Lambda
    out=StabilityCalibNetwork(data=data, pk=pk, Lambda=Lambda, lambda_other_blocks=lambda_other_blocks, pi_list=pi_list, K=K, tau=tau, seed=seed, 
                              method=method, implementation=implementation, resampling_method=resampling_method, scale=scale,
                              PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                              verbose=verbose, debug=FALSE)  
  }
  
  if (verbose){
    cat("\n")
    print("Visited Q:")
    if (nrow(out$Q)>15){
      print(head(out$Q))
      print("[...]")
      print(tail(out$Q))
    } else {
      print(out$Q)
    }
  }
  
  return(out)
}


CalibrateRegression=function(xdata, ydata, Lambda=NULL, pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, 
                             family="gaussian", implementation="glmnet", 
                             resampling_method="subsampling", resampling_function=NULL, PFER_method="MB", PFER_thr=+Inf, FDP_thr=+Inf, 
                             lambda_path_refined_cardinal=100,
                             verbose=TRUE, debug=FALSE, ...){
  if (is.null(Lambda)){
    # Define grid of lambda values (using glmnet implementation)
    Lambda=LambdaGridRegression(xdata=xdata, ydata=ydata, tau=tau, seed=seed, 
                                family=family, implementation=implementation,
                                resampling_method=resampling_method, resampling_function=resampling_function, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                                lambda_path_factor=0.0001, max_density=0.3, 
                                lambda_path_refined_cardinal=lambda_path_refined_cardinal, verbose=verbose, ...)
  }
  
  # Perform stability-based calibration
  out=StabilityCalibRegression(xdata=xdata, ydata=ydata, pk=NULL, Lambda=Lambda, pi_list=pi_list, K=K, tau=tau, seed=seed, 
                               family=family, implementation=implementation,
                               resampling_method=resampling_method, resampling_function=resampling_function, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                               verbose=verbose, debug=debug, ...)
  
  return(out)
}


CalibrationPlot=function(calib_object, block_id=NULL, bi_dim=TRUE, filename=NA, width=10, height=10, lines=TRUE, pt.cex=0.5, lwd=0.2, ymin=NULL, ymax=NULL){
  # Prepare objects
  if ((length(calib_object$params$pk)>1)&(is.null(block_id))){
    bi_dim=FALSE
  }
  if (is.null(ymin)){
    ymin=1
  }
  
  if (bi_dim){
    # Calibration heatmap
    if (is.null(block_id)){
      # ids=which(!is.na(calib_object$Q))
      ids=which(apply(calib_object$S_2d,1,FUN=function(x){any(!is.na(x))}))
      mat=calib_object$S_2d[ids,]
      block_id=1
    } else {
      # ids=which(!is.na(calib_object$Q[,block_id]))
      ids=which(apply(calib_object$S_2d[,,block_id],1,FUN=function(x){any(!is.na(x))}))
      mat=calib_object$S_2d[ids,,block_id]
    }
    rownames(mat)=paste(formatC(calib_object$Lambda[ids,block_id], format="f", digits=4), "-", calib_object$Q[ids,block_id])
    colnames(mat)=calib_object$params$pi_list
    m=matrix("", nrow=nrow(mat), ncol=ncol(mat))
    m[rbind(GetArgmaxId(S=mat))]="X"
    if (is.null(ymax)){
      ymax=nrow(mat)
    }
    mat=mat[ymax:ymin,]
    m=m[ymax:ymin,]
    pheatmap(t(mat), cluster_rows = FALSE, cluster_cols = FALSE, border=NA, display_numbers=t(m),
             filename=filename, width=width, height=height)
  } else {
    # Prepare blocks
    bigblocks=GetBlockMatrix(calib_object$params$pk)
    bigblocks_vect=bigblocks[upper.tri(bigblocks)]
    N_blocks=unname(table(bigblocks_vect))
    blocks=unique(as.vector(bigblocks_vect))
    names(N_blocks)=blocks
    nblocks=max(blocks)
    
    # Calibration plot
    if (!is.na(filename)){
      pdf(filename, width=width, height=height)
    }
    par(mar=c(5,5,3,1), mfrow=c(1,nblocks))
    for (block_id in 1:nblocks){
      # plot(calib_object$Lambda[,block_id], calib_object$S, pch=19, col="navy",
      #      xlab=ifelse(block_id==1, yes=expression(lambda), no=""), 
      #      ylab=ifelse(block_id==1, yes="Stability Score", no=""), 
      #      cex.lab=1.5, cex=pt.cex)
      # abline(v=axTicks(1), lty=2, col="grey")
      # abline(h=axTicks(2), lty=2, col="grey")
      # myS=calib_object$S
      # myS[is.na(myS)]=0
      # points(calib_object$Lambda[,block_id], calib_object$S, pch=19, col="navy", cex=pt.cex)
      # if (lines){
      #   tmp_other_blocks=apply(calib_object$Lambda[,-block_id], 1, paste, collapse="-")
      #   for (i in 1:length(unique(tmp_other_blocks))){
      #     ids=which(tmp_other_blocks==unique(tmp_other_blocks)[i])
      #     if (any(!is.na(calib_object$S[ids]))){
      #       lines(calib_object$Lambda[ids,block_id], calib_object$S[ids], 
      #             col="navy", lwd=0.2, cex=pt.cex)
      #     }
      #   }
      # }
      if (nblocks==1){
        myS_block=calib_object$S
      } else {
        myS_block=calib_object$S_blocks[,block_id]
      }
      plot(calib_object$Lambda[,block_id], myS_block, pch=19, col="navy",
           xlab=ifelse(block_id==1, yes=expression(lambda), no=""),
           ylab=ifelse(block_id==1, yes="Stability Score", no=""),
           cex.lab=1.5, cex=pt.cex)
      abline(v=axTicks(1), lty=2, col="grey")
      abline(h=axTicks(2), lty=2, col="grey")
      myS=myS_block
      myS[is.na(myS)]=0
      points(calib_object$Lambda[,block_id], myS_block, pch=19, col="navy", cex=pt.cex)
      if (lines){
        tmp_other_blocks=apply(calib_object$Lambda[,-block_id], 1, paste, collapse="-")
        for (i in 1:length(unique(tmp_other_blocks))){
          ids=which(tmp_other_blocks==unique(tmp_other_blocks)[i])
          if (any(!is.na(myS_block[ids]))){
            lines(calib_object$Lambda[ids,block_id], myS_block[ids],
                  col="navy", lwd=lwd, cex=pt.cex)
          }
        }
      }
      abline(v=calib_object$Lambda[which.max(myS),block_id], lty=3, col="darkred")
      abline(h=max(myS), lty=3, col="darkred")
    }
    if (!is.na(filename)){
      dev.off()
    }
  }
}


CalibratedAdjacency=function(calib_object, argmax_id=NULL){
  A=matrix(0, ncol=ncol(calib_object$selprop), nrow=nrow(calib_object$selprop))
  bigblocks=GetBlockMatrix(calib_object$params$pk)
  if (is.null(argmax_id)){
    argmax_id=GetArgmaxId(calib_object)
    argmax=GetArgmax(calib_object)
  } else {
    argmax=NULL
    for (block_id in 1:ncol(calib_object$Lambda)){
      argmax=rbind(argmax, c(calib_object$Lambda[argmax_id[block_id,1],block_id],
                             calib_object$params$pi_list[argmax_id[block_id,2]]))
    }
    # print(argmax)
  }
  for (block_id in 1:ncol(calib_object$Lambda)){
    A_block=ifelse(calib_object$selprop[,,argmax_id[block_id,1]]>=argmax[block_id,2],1,0)
    A_block[lower.tri(A_block)]=0
    A_block=A_block+t(A_block) # for symmetry
    if (length(calib_object$params$pk)>1){
      A_block[bigblocks!=block_id]=0
    }
    A=A+A_block
  }
  return(A)
}


CalibratedSelectionProportionsNetwork=function(calib_object, argmax_id=NULL){
  A=matrix(0, ncol=ncol(calib_object$selprop), nrow=nrow(calib_object$selprop))
  bigblocks=GetBlockMatrix(calib_object$params$pk)
  if (is.null(argmax_id)){
    argmax_id=GetArgmaxId(calib_object)
    argmax=GetArgmax(calib_object)
  } else {
    argmax=NULL
    for (block_id in 1:ncol(calib_object$Lambda)){
      argmax=rbind(argmax, c(calib_object$Lambda[argmax_id[block_id,1],],
                             calib_object$params$pi_list[argmax_id[block_id,2]]))
    }
  }
  for (block_id in 1:ncol(calib_object$Lambda)){
    A_block=calib_object$selprop[,,argmax_id[block_id,1]]
    A_block[lower.tri(A_block)]=0
    A_block=A_block+t(A_block) # for symmetry
    if (length(calib_object$params$pk)>1){
      A_block[bigblocks!=block_id]=0
    }
    A=A+A_block
  }
  return(A)
}


CalibratedStableRegression=function(out, argmax_id=NULL){
  if (is.null(argmax_id)){
    argmax_id=GetArgmaxId(out)
  } 
  stability_selected=ifelse(out$selprop[argmax_id[1],]>=out$params$pi_list[argmax_id[2]],1,0)
  return(stability_selected)
}


CalibratedSelectionProportionsRegression=function(out, filename=NULL, width=15, height=8, lambda=NULL, annot=NULL, plot=TRUE){
  # stability_selected=Calib(out)
  # names(stability_selected[stability_selected==1])
  
  if (!is.null(filename)){
    pdf(filename, width=width, height=height)
  }
  par(mar=c(1,5,1,1))
  # if (is.null(lambda)){
  hat_theta=GetArgmaxId(out)
  m=out$selprop[hat_theta[1],]
  calibrated_pi=out$params$pi_list[hat_theta[2]]
  # } else {
  #   m=out$stab[lambda,]
  #   calibrated_pi=10
  # }
  m=m[sort.list(m, decreasing=TRUE)]
  if (!is.null(annot)){
    names(m)=annot[names(m)]
  }
  if (plot){
    mycolours=ifelse(m>=(1-calibrated_pi), yes=1, no=0)+ifelse(m>=calibrated_pi, yes=1, no=0)
    plot(m, type="h", lwd=1, xaxt="n", las=1,ylim=c(0,1),
         xlab="", ylab="Selection Proportion", cex.lab=1.5,
         col=c("darkgrey", "red", "darkred")[mycolours+1])
    abline(h=calibrated_pi, col="darkred", lty=2)
    abline(h=1-calibrated_pi, col="darkred", lty=2)
  }
  # axis(side=1, at=1:length(m), labels=names(m), las=2)
  if (!is.null(filename)){
    dev.off()
  }
  return(m)
}


GetGraph=function(calib_object=NULL, adjacency=NULL, 
                  node_label=NULL, node_color=NULL, node_shape=NULL,
                  weighted=NULL, satellites=FALSE){
  # either out or adjacency have to be provided
  
  if (is.null(adjacency)){
    if (is.null(calib_object)){
      stop("Either 'calib_object' or 'adjacency' needs to be provided.")
    }
    adjacency=CalibratedAdjacency(calib_object)
  }
  
  if (is.null(node_color)){
    node_color=rep("skyblue", ncol(adjacency))
  }
  
  if (is.null(node_shape)){
    node_shape=rep("circle", ncol(adjacency))
  }
  
  if (is.null(node_label)){
    node_label=colnames(adjacency)
  }
  
  names(node_color)=colnames(adjacency)
  names(node_label)=colnames(adjacency)
  names(node_shape)=colnames(adjacency)
  
  mygraph=graph_from_adjacency_matrix(adjacency, mode = 'undirected', weighted = weighted)
  
  if (!satellites){
    mygraph=delete.vertices(mygraph, v=names(degree(mygraph))[degree(mygraph)==0])
  }
  
  V(mygraph)$size=as.numeric(as.character(cut(degree(mygraph), breaks = 4, labels = c(3, 4, 5, 6))))
  V(mygraph)$label=node_label[V(mygraph)$name]
  V(mygraph)$color=node_color[V(mygraph)$name]
  V(mygraph)$shape=node_shape[V(mygraph)$name]
  V(mygraph)$frame.color=V(mygraph)$color
  V(mygraph)$label.family="sans"
  E(mygraph)$color="grey60"
  V(mygraph)$label.cex=as.numeric(as.character(cut(degree(mygraph), breaks = 4, labels = c(0.4, 0.45, 0.5, 0.55))))
  V(mygraph)$label.color="grey20"
  E(mygraph)$width=0.5
  
  if (!is.null(weighted)){
    E(mygraph)$color=c('red', 'blue', 'forestgreen')[E(mygraph)$weight]
  }
  
  return(mygraph)
}


GetPerformance=function(A, cor=NULL, thr=0.5){
  
  if (is.matrix(A)){
    p=ncol(A)
    N=p*(p-1)/2
    A=A[upper.tri(A)]
  } else {
    N=length(A)
  }
  
  TP=sum(A==3)
  FN=sum(A==2)
  FP=sum(A==1)
  TN=sum(A==0)
  
  if (!is.null(cor)){
    if (is.matrix(cor)){
      cor_vect=cor[upper.tri(cor)]
    } else {
      cor_vect=cor
    }
    FP_c=sum((A==1)&(abs(cor_vect)>=thr))
    FP_i=sum((A==1)&(abs(cor_vect)<thr))
  }
  
  sensitivity=TP/(TP+FN)
  specificity=TN/(TN+FP)
  accuracy=(TP+TN)/N
  if (TP+FP>0){
    precision=TP/(TP+FP)
  } else {
    precision=0
  }
  if ((TP+FN)>0){
    recall=TP/(TP+FN)
  } else {
    recall=1
  }
  if ((precision>0)|(recall>0)){
    F1_score=2*precision*recall/(precision+recall)
  } else {
    F1_score=0
  }
  
  if (is.null(cor)){
    return(data.frame(TP=TP, FN=FN, FP=FP, TN=TN, 
                      sensitivity=sensitivity, specificity=specificity, 
                      accuracy=accuracy, precision=precision, recall=recall, F1_score=F1_score))
  } else {
    return(data.frame(TP=TP, FN=FN, FP=FP, TN=TN, FP_c=FP_c, FP_i=FP_i,
                      sensitivity=sensitivity, specificity=specificity, 
                      accuracy=accuracy, precision=precision, recall=recall, F1_score=F1_score))
  }
}


PerformanceGraphRepresentation=function(Asum, node_color=NULL, node_shape=NULL, plot=FALSE, filename=NULL, width=7, height=7, 
                                        mycolours=c("grey95", "tomato", "forestgreen", "navy"), mylty=c(4,2,3,1)){
  # Refine inputs
  names(mycolours)=names(mylty)=c("TN", "FP", "FN", "TP")
  if (is.null(node_color)){
    node_color=rep("black",ncol(Asum))
  }
  
  # Make consensus graph
  g=GetGraph(adjacency=ifelse(Asum!=0,yes=1,no=0), node_color=node_color, node_shape=node_shape)
  
  # Format edges
  V(g)$size=V(g)$size/3+1
  V(g)$label=rep("", length(V(g)$label))
  myedgecolour=mycolours[Asum[get.edgelist(g)]+1]
  myedgelty=c(2, 3, 1)[Asum[get.edgelist(g)]]
  E(g)$color=myedgecolour
  E(g)$width=1
  E(g)$lty=myedgelty
  
  # Plot graph
  if (plot|(!is.null(filename))){
    if (!is.null(filename)){
      pdf(filename, width=width, height=height)
    }
    par(mar=rep(0,4))
    plot(g, layout=layout_with_kk(g))
    if (!is.null(filename)){
      dev.off()
    }
  }
  
  # Return output graph
  return(g)
}


GetBlockPerformanceNetwork=function(Asum, pk, cor=NULL){
  Asum_vect=Asum[upper.tri(Asum)]
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  if (!is.null(cor)){
    cor_vect=cor[upper.tri(cor)]
  } else {
    cor_vect=NULL
  }
  
  return(rbind(GetPerformance(Asum, cor=cor),
               GetPerformance(Asum_vect[bigblocks_vect==1], cor=cor_vect[bigblocks_vect==1]),
               GetPerformance(Asum_vect[bigblocks_vect==2], cor=cor_vect[bigblocks_vect==2]),
               GetPerformance(Asum_vect[bigblocks_vect==3], cor=cor_vect[bigblocks_vect==3])))
}


MergeOutputsNetwork=function(output1, output2, ordered=TRUE){
  # Checking consistency between both outputs
  if(!all(unlist(output1$methods)==unlist(output2$methods))){
    print("output1:")
    print(unlist(output1$methods))
    print("output2:")
    print(unlist(output2$methods))
    stop(paste0("Different methods were used for these outputs, they cannot be merged."))
  }
  
  if(!all(unlist(output1$params)==unlist(output2$params))){
    print("output1:")
    print(unlist(output1$params))
    print("output2:")
    print(unlist(output2$params))
    stop(paste0("Different parameters were used for these outputs, they cannot be merged."))
  }
  
  if ((nrow(output1$selprop)!=nrow(output2$selprop))|(ncol(output1$selprop)!=ncol(output2$selprop))){
    stop("Different numbers of nodes were used for these outputs, they cannot be merged.")
  }
  
  # Merging the entries
  out=NULL
  out$S=rbind(output1$S, output2$S)
  out$Lambda=rbind(output1$Lambda, output2$Lambda)
  out$Q=rbind(output1$Q, output2$Q)
  out$Q_s=rbind(output1$Q_s, output2$Q_s)
  out$P=rbind(output1$P, output2$P)
  out$PFER=rbind(output1$PFER, output2$PFER)
  out$FDP=rbind(output1$FDP, output2$FDP)
  out$S_2d=rbind(output1$S_2d, output2$S_2d)
  out$PFER_2d=rbind(output1$PFER_2d, output2$PFER_2d)
  out$FDP_2d=rbind(output1$FDP_2d, output2$FDP_2d)
  
  # Merging the selection proportions
  out$selprop=abind(output1$selprop, output2$selprop, along=3)
  
  # Filling the methods and params
  out$methods=output1$methods
  out$params=output1$params
  
  # Order the output if required
  if (ordered){
    out=OrderOutputNetwork(out)
  }
  
  # Return the output
  return(out)
}


OrderOutputNetwork=function(output){
  if (ncol(output$Q)>1){
    stop("This function cannot be used to order outputs of the multi-block calibration.")
  }
  
  # Ordering the entries based on Lambda
  Lambda=output$Lambda[,1]
  output$S=output$S[sort.list(Lambda, decreasing=TRUE),,drop=FALSE]
  output$Q=output$Q[sort.list(Lambda, decreasing=TRUE),,drop=FALSE]
  output$Q_s=output$Q_s[sort.list(Lambda, decreasing=TRUE),,drop=FALSE]
  output$P=output$P[sort.list(Lambda, decreasing=TRUE),,drop=FALSE]
  output$PFER=output$PFER[sort.list(Lambda, decreasing=TRUE),,drop=FALSE]
  output$FDP=output$FDP[sort.list(Lambda, decreasing=TRUE),,drop=FALSE]
  output$S_2d=output$S_2d[sort.list(Lambda, decreasing=TRUE),,drop=FALSE]
  output$PFER_2d=output$PFER_2d[sort.list(Lambda, decreasing=TRUE),,drop=FALSE]
  output$FDP_2d=output$FDP_2d[sort.list(Lambda, decreasing=TRUE),,drop=FALSE]
  
  # Ordering the selection proportions
  output$selprop=output$selprop[,,sort.list(Lambda, decreasing=TRUE)]
  
  # Ordering Lambda
  output$Lambda=output$Lambda[sort.list(Lambda, decreasing=TRUE),,drop=FALSE]
  
  return(output)
}


SelectionFunction=function(x, y, lambda, family, implementation="glmnet", ...){
  # implementation: name of a function (character string) with a variable/group selection algorithm which outputs a matrix 
  # of coefficients (beta/loadings or simply binary indicators) where columns are variables and rows are parameter values
  
  # Making sure none of the variables has a null standard deviation
  mysd=apply(x,2,sd)
  if (any(mysd==0)){
      for (k in which(mysd==0)){
          x[,k]=x[,k]+rnorm(n=nrow(x), sd=mysd[k]/100)
      }
  }
  x=scale(x)
  
  if (implementation=="glmnet"){
    # Running the regression
    mymodel=glmnet(x=x, y=y, lambda=lambda, family=family, ...)
    
    # Extracting and formatting the beta coefficients
    if (family!="mgaussian"){
      mybeta=coef(mymodel)
      mybeta=t(as.matrix(mybeta))
      mybeta=mybeta[,colnames(x)] # removing the intercept if included
    } else {
      mybeta=array(NA, dim=c(length(lambda), ncol(x), ncol(y)), 
                   dimnames=list(paste0("s",0:(length(lambda)-1)), colnames(x), colnames(y)))
      for (y_id in 1:ncol(y)){
        tmpbeta=coef(mymodel)[[y_id]]
        tmpbeta=t(as.matrix(tmpbeta))
        tmpbeta=tmpbeta[,colnames(x),drop=FALSE] # removing the intercept if included
        mybeta[rownames(tmpbeta),colnames(tmpbeta),y_id]=tmpbeta
      }
    }
  }
  
  if (implementation=="gglasso"){
    # Running the regression
    if (family=="binomial"){
      ytmp=y
      ytmp[ytmp==min(ytmp)]=-1
      ytmp[ytmp==max(ytmp)]=1
      mymodel=gglasso(x, ytmp, lambda=lambda, loss="logit", ...)
    }
    if (family=="gaussian"){
      mymodel=gglasso(x, y, lambda=lambda, loss="ls", ...)
    }
    if (!family%in%c("binomial", "gaussian")){
      stop("Please set family as 'binomial' or 'gaussian' for use with gglasso.")
    }
    
    # Extracting and formatting the beta coefficients
    mybeta=t(as.matrix(mymodel$beta))
    mybeta=mybeta[,colnames(x)]
  }
  
  if (!implementation%in%c("glmnet", "gglasso")){
    mybeta=do.call(get(implementation), args=list(x=x, y=y, lambda=lambda, family=family, ...))
  }
  
  return(mybeta)
}


huge.adjacency.generator=function(d = 50, graph = "random", g = NULL, prob = NULL, verbose = FALSE){
  if (verbose) 
    cat("Generating data from the multivariate normal distribution with the", 
        graph, "graph structure....")
  if (is.null(g)) {
    g = 1
    if (graph == "hub" || graph == "cluster") {
      if (d > 40) 
        g = ceiling(d/20)
      if (d <= 40) 
        g = 2
    }
  }
  if (graph == "random") {
    if (is.null(prob)) 
      prob = min(1, 3/d)
    prob = sqrt(prob/2) * (prob < 0.5) + (1 - sqrt(0.5 - 0.5 * prob)) * (prob >= 0.5)
  }
  if (graph == "cluster") {
    if (is.null(prob)) {
      if (d/g > 30) 
        prob = 0.3
      if (d/g <= 30) 
        prob = min(1, 6 * g/d)
    }
    prob = sqrt(prob/2) * (prob < 0.5) + (1 - sqrt(0.5 - 
                                                     0.5 * prob)) * (prob >= 0.5)
  }
  if (graph%in%c("cluster", "hub")){
    g.large = d%%g
    g.small = g - g.large
    n.small = floor(d/g)
    n.large = n.small + 1
    g.list = c(rep(n.small, g.small), rep(n.large, g.large))
    g.ind = rep(c(1:g), g.list)
    rm(g.large, g.small, n.small, n.large, g.list)
    gc()
  }
  theta = matrix(0, d, d)
  if (graph == "band") {
    for (i in 1:g) {
      diag(theta[1:(d - i), (1 + i):d]) = 1
      diag(theta[(1 + i):d, 1:(d - 1)]) = 1
    }
  }
  if (graph == "cluster") {
    for (i in 1:g) {
      tmp = which(g.ind == i)
      tmp2 = matrix(runif(length(tmp)^2, 0, 0.5), length(tmp), 
                    length(tmp))
      tmp2 = tmp2 + t(tmp2)
      theta[tmp, tmp][tmp2 < prob] = 1
      rm(tmp, tmp2)
      gc()
    }
  }
  if (graph == "hub") {
    for (i in 1:g) {
      tmp = which(g.ind == i)
      theta[tmp[1], tmp] = 1
      theta[tmp, tmp[1]] = 1
      rm(tmp)
      gc()
    }
  }
  if (graph == "random") {
    tmp = matrix(runif(d^2, 0, 0.5), d, d)
    tmp = tmp + t(tmp)
    theta[tmp < prob] = 1
    rm(tmp)
    gc()
  }
  if (graph == "scale-free") {
    out = .Call("_huge_SFGen", 2, d)
    theta = matrix(as.numeric(out$G), d, d)
    ids=sample(ncol(theta))
    theta=theta[ids,ids] # re-organising to avoid decreasing proba over higher variable index
  }
  diag(theta) = 0
  colnames(theta)=rownames(theta)=paste0("var",1:ncol(theta))
  
  return(theta)
}


SimulateGraph=function(n=100, pk=10, u=NULL, v_within=1, v_between=0.5, 
                       topology="random", nu=0.1, output_matrices=FALSE, pd_strategy="diagonally_dominance",
                       all_negative=FALSE){
  # pd_strategy: "diagonally_dominance" or "nonnegative_eigenvalues"
  
  # Define grid of u values if not provided
  if (is.null(u)){
    u=10^-(seq(0,5,by=0.1))
    refining_u_grid=TRUE
    niter_max=5
  } else {
    refining_u_grid=FALSE
  }
  
  # Define number of nodes
  p=sum(pk)
  
  # Create matrix with block indices
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  N_blocks=unname(table(bigblocks_vect))
  block_ids=unique(as.vector(bigblocks))
  names(N_blocks)=block_ids
  nblocks=max(block_ids)
  
  # Build v matrix
  v_list=rep(NA,length(block_ids))
  v_list[diag(bigblocks)]=v_within
  v_list[is.na(v_list)]=v_between
  v=bigblocks
  for (k in block_ids){
    v[bigblocks==k]=v_list[k]
  }
  
  # Simulation of the adjacency matrix
  if ((topology=="random")&(nblocks>1)){
    # # Ensuring that the density is homogeneous across blocks for multi-block simulation
    # theta_vect=rep(0, length(bigblocks_vect))
    # for (block_id in 1:nblocks){
    #   tmp=rbinom(n=N_blocks[block_id], size=1, prob=nu)
    #   theta_vect[bigblocks_vect==block_id]=tmp
    # }
    # theta=matrix(0, ncol=p, nrow=p)
    # theta[upper.tri(theta)]=theta_vect
    # theta=theta+t(theta)
    # rownames(theta)=colnames(theta)=paste0("var", 1:p)
    
    # Formula from huge package
    nu = sqrt(nu/2) * (nu < 0.5) + (1 - sqrt(0.5 - 0.5 * nu)) * (nu >= 0.5)
    # Ensuring that the density is homogeneous across blocks for multi-block simulation
    theta_vect=rep(0, length(bigblocks_vect))
    for (block_id in 1:nblocks){
      tmp=runif(N_blocks[block_id], 0, 0.5)+runif(N_blocks[block_id], 0, 0.5)
      tmp=ifelse(tmp<nu, yes=1, no=0)
      theta_vect[bigblocks_vect==block_id]=tmp
    }
    theta=matrix(0, ncol=p, nrow=p)
    theta[upper.tri(theta)]=theta_vect
    theta=theta+t(theta)
    rownames(theta)=colnames(theta)=paste0("var", 1:p)
  } else {
    # Using the huge function for other topologies
    theta=huge.adjacency.generator(d=p, prob=nu, graph=topology)
  }
  
  # Fill off-diagonal entries of the precision matrix
  omega = theta * v
  # omega_within = theta * v_within
  
  # Calibrate u based on contrasts of the correlation matrix
  contrast=NULL
  for (u_value in u){
    omega_tmp=omega
    if (pd_strategy=="diagonally_dominance"){
      diag(omega_tmp) = apply(abs(omega_tmp),1,sum) + u_value
    }
    if (pd_strategy=="nonnegative_eigenvalues"){
      diag(omega_tmp) = abs(min(eigen(omega_tmp)$values)) + u_value # used in huge
    }
    C = cov2cor(solve(omega_tmp))
    contrast=c(contrast, GetContrast(C))
  }
  
  # Avoid extreme values in u grid if not provided by the user
  tolerance=10
  if (refining_u_grid){
    stop=0
    niter=1
    while (stop==0){
      niter=niter+1
      if (niter==niter_max){
        stop=1
      }
      if (any(which(contrast==max(contrast))%in%seq(tolerance,length(u)-tolerance)==TRUE)){
        stop=1
      } else {
        if (any(which(contrast==max(contrast))%in%seq(1,tolerance)==TRUE)){
          u=c(u, 10^-seq(min(-log10(u))-5, min(-log10(u)), by=0.1))
        }
        if (any(which(contrast==max(contrast))%in%seq(length(u)-tolerance,length(u))==TRUE)){
          u=c(u, 10^-seq(max(-log10(u)), max(-log10(u)+5), by=0.1))
        }
        u=sort(u, decreasing=TRUE)
        contrast=NULL
        for (u_value in u){
          omega_tmp=omega
          if (pd_strategy=="diagonally_dominance"){
            diag(omega_tmp) = apply(abs(omega_tmp),1,sum) + u_value
          }
          if (pd_strategy=="nonnegative_eigenvalues"){
            diag(omega_tmp) = abs(min(eigen(omega_tmp)$values)) + u_value # used in huge
          }
          C = cov2cor(solve(omega_tmp))
          contrast=c(contrast, GetContrast(C))
        }
      }
    }
  }
  
  # Compute calibrated precision matrix
  if (length(u)>1){
    u_value=u[which.max(contrast)]
    if (pd_strategy=="diagonally_dominance"){
      diag(omega) = apply(abs(omega),1,sum) + u_value
    }
    if (pd_strategy=="nonnegative_eigenvalues"){
      diag(omega) = abs(min(eigen(omega)$values)) + u_value # used in huge
    }
  } else {
    omega=omega_tmp
  }
  
  if (!all_negative){
    # Simulate the sign of the associations (uniform)
    sign_mat=matrix(0, nrow=nrow(omega), ncol=ncol(omega))
    sign_mat[upper.tri(sign_mat)]=sample(c(-1,1), size=sum(upper.tri(omega)), replace=TRUE)
    sign_mat=sign_mat+t(sign_mat)
    diag(sign_mat)=1
    omega=omega*sign_mat
  }
  
  # Compute the correlation matrix
  C = cov2cor(solve(omega)) # true correlation matrix - called sigma in huge
  
  # Compute the partial correlation matrix
  if (output_matrices){
    phi=-cov2cor(omega)+2*diag(ncol(omega))
  }
  
  # Simulate data from multivariate normal distribution
  x = mvrnorm(n, rep(0, p), C)
  colnames(x)=paste0("var",1:ncol(x))
  rownames(x)=paste0("obs",1:nrow(x))
  
  if (output_matrices){
    return(list(data=x, theta=theta, 
                omega=omega, phi=phi, C=C,
                u=u_value, u_grid=u, contrast_path=contrast))
  } else {
    return(list(data=x, theta=theta))
  }
}


SimulateXY=function(n=100, pk=10, X=NULL, underlying_graph=FALSE, u=10^-(seq(0,5,by=0.1)), v_within=1, v_between=0.5, 
                    topology="random", nu=0.1, nu_pred=0.2, beta_min=1, beta_max=1, sd_pred_error=1, family="gaussian",
                    output_matrices=FALSE, pd_strategy="diagonally_dominance", all_negative=FALSE){
  p=sum(pk)
  
  if (underlying_graph){ # DO NOT USE, UNDER DEVELOPMENT
    # Simulate the predictor matrix X
    simul=SimulateGraph(n=n, pk=pk, u=u, v_within=v_within, v_between=v_between,
                        topology=topology, nu=nu, output_matrices=output_matrices, pd_strategy=pd_strategy,
                        all_negative=all_negative)
    X=simul$data
    
    # Select the hub
    g=GetGraph(adjacency=simul$theta)
    theta_pred=rep(0,ncol(simul$data))
    names(theta_pred)=colnames(simul$data)
    theta_pred[which.max(degree(g))]=1
    
    # Select the neighbours with probability exponentially decreasing with the shortest path
    mysp=shortest.paths(g)
    myprob=exp(max(mysp[which(theta_pred==1),])-mysp[which(theta_pred==1),])
    myprob=(myprob)/sum(myprob)
    theta_pred[sample(1:ncol(simul$data), size=ceiling(nu_pred*ncol(simul$data)), prob=myprob)]=1
    true_predictors=which(theta_pred==1)
    
    p=ncol(X)
    n=nrow(X)
  } 
  
  if ((is.null(X))&(!underlying_graph)) {
    X=NULL
    for (k in 1:p){
      X=cbind(X, rnorm(n, mean=0, sd=1))
    }
  }
  
  colnames(X)=paste0("var",1:ncol(X))
  rownames(X)=paste0("obs",1:nrow(X))
  
  # Get the binary vector of true predictors
  theta_pred=rbinom(p, size=1, prob=nu_pred)
  names(theta_pred)=colnames(X)
  
  # Simulate a vector of betas
  # beta=runif(p, min=beta_min, max=beta_max)*sample(c(-1,1), size=p, replace=TRUE)
  beta=runif(p, min=beta_min, max=beta_max)
  beta=beta*theta_pred
  
  # Compute the predicted values of Y
  Y_pred=X%*%beta
  
  # Introduce some gaussian error
  Y=Y_pred+rnorm(n, mean=0, sd=sd_pred_error)
  
  # Compute binary outcome for logistic regression
  if (family=="binomial"){
    proba=1/(1+exp(-Y)) # inverse logit 
    Y_bin=rbinom(n, size=1, prob=proba)
  }
  
  # Return the simulated X and Y
  if (family=="binomial"){
    out=list(X=X, Y=Y_bin, proba=proba, logit_proba=Y, logit_proba_pred=Y_pred, theta_pred=theta_pred)
  } else {
    out=list(X=X, Y=Y, Y_pred=Y_pred, theta_pred=theta_pred)
  }
  if (underlying_graph){
    simul$data=NULL
    return(c(out, simul))
  } else {
    return(out)
  }
}


GetBlockDistributions=function(mat, pk, filename=NULL){
  mat_vect=mat[upper.tri(mat)]
  
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  blocks=sort(unique(bigblocks_vect))
  
  for (k in blocks){
    tmpmat=mat_vect[bigblocks_vect==k]
    tmpd=density(tmpmat, from=-1, to=1)
    tmpsd=sd(tmpmat)
    assign(paste0("mat",k), tmpmat)
    assign(paste0("d",k), tmpd)
    assign(paste0("sd",k), tmpsd)
  }
  
  mycolors=MyPalPairs=brewer.pal(n=min(12,length(blocks)),name='Paired')
  mycolors=colorRampPalette(mycolors)(length(blocks))
  
  if (!is.null(filename)){
    pdf(filename)
  }
  par(mar=c(5,1,1,1))
  plot(d1, col=mycolors[1], ylim=range(c(d1$y, d2$y, d3$y)), lwd=2,
       panel.first=abline(v=0, lty=2, col="black"), las=1, 
       xlab="Pearson's correlation", main="", yaxt="n", ylab="", cex.lab=1.5)
  for (k in blocks[2:length(blocks)]){
    tmpd=eval(parse(text=paste0("d",k)))
    lines(tmpd, col=mycolors[k], lwd=2)
  }
  legend("topleft", lty=1, lwd=2, col=mycolors,
         legend=sapply(blocks, FUN=function(k){paste0("Block ",k, ": ", 
                                                      ifelse(k%in%unique(diag(bigblocks)), yes="within", no="between"),
                                                      " (sd=", round(eval(parse(text=paste0("sd",k))), digits=3), ")")})
  )
  if (!is.null(filename)){
    dev.off()
  }
  
  output=sapply(blocks, FUN=function(k){eval(parse(text=paste0("sd",k)))})
  names(output)=paste0("sd", blocks)
  
  return(output)
}


GetPRDistance=function(precision, recall){
  # dist_visited=1/2*(2-((1-precision)^2+(1-recall)^2))
  dist_visited=sqrt((1-precision)^2+(1-recall)^2) # Euclidian distance to (1,1) point
  return(list(distance=dist_visited, rank=floor(rank(dist_visited))))
}


glasso.graphical_model=function (x, y, q, scale=TRUE, ...){
  if (!requireNamespace("QUIC")) 
    stop("Package ", sQuote("QUIC"), " is required but not available")
  extraargs <- list(...)
  if (scale){
    empirical.cov <- cor(x)
  } else {
    empirical.cov <- cov(x)
  }
  lams <- extraargs$lams
  if (is.null(lams)) {
    max.cov <- max(abs(empirical.cov[upper.tri(empirical.cov)]))
    lams <- getLamPath(max.cov, max.cov * 0.05, len = 40)
  }
  est=NULL
  est$X=NULL
  for (k in 1:length(lams)){
    # est <- QUIC::QUIC(empirical.cov, rho = 1, path = lams, msg = 0)
    myglasso=glassoFast(empirical.cov, rho=lams[k])
    est$X=abind(est$X, myglasso$wi, along=3)
  }
  ut <- upper.tri(empirical.cov)
  qvals <- sapply(1:length(lams), function(idx) {
    m <- est$X[, , idx]
    sum(m[ut] != 0)
  })
  qq <- qvals >= q
  if (!any(qq)) {
    stop("Didn't reach the required number of variables. Try supplying lambda manually")
  }
  lamidx <- which.max(qvals >= q)
  M <- est$X[, , lamidx][ut]
  selected <- (M != 0)
  s <- sapply(1:lamidx, function(idx) {
    m <- est$X[, , idx][ut] != 0
    return(m)
  })
  colnames(s) <- as.character(1:ncol(s))
  return(list(selected = selected, path = s))
}

class(glasso.graphical_model)=c("function", "graphical_model")


glasso.pulsar=function (data, lambda, scale=TRUE){
  x=data
  # extraargs <- list(...)
  
  if (scale){
    empirical.cov <- cor(x)
  } else {
    empirical.cov <- cov(x)
  }
  lams <- lambda
  if (is.null(lams)) {
    max.cov <- max(abs(empirical.cov[upper.tri(empirical.cov)]))
    lams <- getLamPath(max.cov, max.cov * 0.05, len = 40)
  }
  path=list()
  for (k in 1:length(lams)){
    # est <- QUIC::QUIC(empirical.cov, rho = 1, path = lams, msg = 0)
    myglasso=glassoFast(empirical.cov, rho=lams[k])
    A=ifelse(myglasso$wi!=0, yes=1, no=0)
    A=A+t(A)
    A=ifelse(A!=0, yes=1, no=0)
    path=c(path, list(A))
  }
  return(list(path=path))
}

glmnet.lasso_model=function(x, y, q, type = c("conservative", "anticonservative"), ...){
  if (!requireNamespace("glmnet", quietly = TRUE)) 
    stop("Package ", sQuote("glmnet"), " needed but not available")
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"), " is coerced to a model matrix without intercept")
    x <- model.matrix(~. - 1, x)
  }
  type <- match.arg(type)
  if (type == "conservative") 
    fit <- suppressWarnings(glmnet::glmnet(x, y, pmax = q, 
                                           ...))
  if (type == "anticonservative") 
    fit <- glmnet::glmnet(x, y, dfmax = q - 1, ...)
  selected <- predict(fit, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  cf <- fit$beta
  sequence <- as.matrix(cf != 0)
  return(list(selected = ret, path = sequence))
}



