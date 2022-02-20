

### Simulating data
sim_train <- function(pre, delta, std, dvname, logrt=FALSE, N_subj=30, rho=0.5, seed=2022){
  
  # pre: mean of pre-test performance (assuming the same for the two groups);
  # delta: the change from pre- to post-tests;
  # std: the same standard deviation were used for all conditions;
  # dvname: the name of the dependent variable;
  # N_subj: number of participant in each group;
  # rho: correlation between pre- and post-tests.
  
  set.seed(seed)
  
  Groups <- c("control", "train") # Between-subject factor
  Test <- c("post_test", "pre_test") # Within-subject factor
  
  # construct the variance-covariance matrix
  C <- matrix(rho, nrow = 2, ncol = 2)
  diag(C) <- 1
  
  post <- delta + pre
  
  control <- mvtnorm::rmvnorm(N_subj, c(pre, post[1]), sigma = C*std^2)
  train <- mvtnorm::rmvnorm(N_subj, c(pre, post[2]), sigma = C*std^2)
  
  df_sim <- rbind(control, train) %>% 
    as_tibble() %>% 
    transmute(Subj = 1:(N_subj*length(Groups)),
              Group = rep(Groups, each = N_subj),
              pre = V1,
              post = V2) %>% 
    pivot_longer(c(pre, post), names_to = "Test", values_to = dvname)
  
  if (logrt) {
    df_sim <- df_sim %>% 
      mutate(rt = exp(.data[[dvname]]))
  }
  
  return(df_sim)
}


### Simulating HBER for the hypothetical training study
sim_hber_null <- function (N_subj = 30, N_iter = 1000, thealpha = 0.05, 
                           seed = 2022, file_cache = NULL) {
  
  if (!is.null(file_cache) && file.exists(as.character(file_cache))){
    out <- read_rds(file_cache)
    
  } else {
    
    p_d_simp <- rep(NaN, N_iter)
    sig_d_simp <- rep(NaN, N_iter)
    p_d_inter <- rep(NaN, N_iter)
    sig_d_inter <- rep(NaN, N_iter)
    
    p_rt_simp <- rep(NaN, N_iter)
    sig_rt_simp <- rep(NaN, N_iter)
    p_rt_inter <- rep(NaN, N_iter)
    sig_rt_inter <- rep(NaN, N_iter)
    
    conclusion <- rep(NaN, N_iter)
    
    # Each simulation
    for (i in 1:N_iter) {
      # simulate null data for d
      df_d <- sim_train(pre=2.8, delta=c(0,0), std=0.5, dvname="d", 
                        logrt=F, N_subj=N_subj, rho=0.5, seed=seed+i)
      
      # simulate null data for logrt (rt)
      df_rt <- sim_train(pre=6.4, delta=c(0,0), std=0.1, dvname="logrt", 
                         logrt=T, N_subj=N_subj, rho=0.5, seed=seed+i)
      
      # ANOVA for d
      d_aov <- aov_4(d ~ Group * Test + (Test | Subj), df_d)
      d_emm <- emmeans(d_aov, ~ Group + Test)
      
      df_d_simp <- as_tibble(contrast(d_emm, "revpairwise", by="Group")[2])
      p_d_simp[i] <- df_d_simp[["p.value"]]
      sig_d_simp[i] <- p_d_simp[i] < thealpha
      
      df_d_inter <- as_tibble(contrast(d_emm, interaction="revpairwise"))
      p_d_inter[i] <- df_d_inter[["p.value"]]
      sig_d_inter[i] <- p_d_inter[i] < thealpha
      
      # ANOVA for logrt
      rt_aov <- aov_4(logrt ~ Group * Test + (Test | Subj), df_rt)
      rt_emm <- emmeans(rt_aov, ~ Group + Test)
      
      df_rt_simp <- as_tibble(contrast(rt_emm, "revpairwise", by="Group")[2])
      p_rt_simp[i] <- df_rt_simp[["p.value"]]
      sig_rt_simp[i] <- p_rt_simp[i] < thealpha
      
      df_rt_inter <- as_tibble(contrast(rt_emm, interaction="revpairwise"))
      p_rt_inter[i] <- df_rt_inter[["p.value"]]
      sig_rt_inter[i] <- p_rt_inter[i] < thealpha
      
      # EOI
      conclusion[i] <- (sig_d_simp[i] & sig_d_inter[i]) | 
        (sig_rt_simp[i] & sig_rt_inter[i])
    }
    
    out <- tibble(iteration = 1:N_iter, conclusion,
                  p_d_simp, sig_d_simp, p_d_inter, sig_d_inter,
                  p_rt_simp, sig_rt_simp, p_rt_inter, sig_rt_inter)
    
    if (!is.null(file_cache)){
      write_rds(out, file=file_cache)
    }
  }
  
  return(out)
}
