

### Simulating data
sim_train <- function(pre, delta, std, dvname, logrt=FALSE, N_subj=30, rho=0.5, seed=2022){
  
  # pre: mean of pre-test performance (assuming the same for the two groups);
  # delta: the change from pre- to post-tests (control, train);
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
  
  if (length(pre)==1){
    pre = rep(pre, 1, 2)
  }
  
  post <- delta + pre
  
  control <- mvtnorm::rmvnorm(N_subj, c(pre[1], post[1]), sigma = C*std^2)
  train <- mvtnorm::rmvnorm(N_subj, c(pre[2], post[2]), sigma = C*std^2)
  
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
sim_train_null <- function (iter = 1000, N_subj = 30, N_core = 2, 
                           file_cache = NULL) {
  # This function simulates data, performs ANOVA, and records the p-values for 
  # one simple effect (post-test vs. pre-test in the training group) and the 
  # two-way interaction (performance increase in training group relative to the
  # control group).
  
  # iter: iteration number for each sample size
  # N_subj: sample size for each group.
  # N_core: number of cores to run the simulation
  # file_cache: file name of the cache file. If NULL, no file will be cached.
  
  # Simulation for one iteration
  sim_hber_single <- function(N_s) {
    # simulate null data for d
    df_d <- sim_train(pre=2.8, delta=c(0,0), std=0.5, dvname="d", 
                      logrt=F, N_subj=N_s, rho=0.5, seed=NULL)
    
    # simulate null data for logrt (rt)
    df_rt <- sim_train(pre=6.4, delta=c(0,0), std=0.1, dvname="logrt", 
                       logrt=T, N_subj=N_s, rho=0.5, seed=NULL)
    
    # ANOVA for d
    d_aov <- aov_4(d ~ Group * Test + (Test | Subj), df_d)
    d_emm <- emmeans(d_aov, ~ Group + Test)
    
    df_d_simp <- as_tibble(contrast(d_emm, "revpairwise", by="Group")[2])
    p_d_simp <- df_d_simp[["p.value"]] # the simple effect
    df_d_inter <- as_tibble(contrast(d_emm, interaction="revpairwise"))
    p_d_inter <- df_d_inter[["p.value"]] # the interaction
    
    # ANOVA for logrt
    rt_aov <- aov_4(logrt ~ Group * Test + (Test | Subj), df_rt)
    rt_emm <- emmeans(rt_aov, ~ Group + Test)
    
    df_rt_simp <- as_tibble(contrast(rt_emm, "revpairwise", by="Group")[2])
    p_rt_simp <- df_rt_simp[["p.value"]] # the simple effect
    df_rt_inter <- as_tibble(contrast(rt_emm, interaction="revpairwise"))
    p_rt_inter <- df_rt_inter[["p.value"]] # the interaction
    
    rawp <- tibble(N_subj = N_s,
                   p_d_simple = p_d_simp,
                   p_d_inter = p_d_inter,
                   p_rt_simple = p_rt_simp,
                   p_rt_inter = p_rt_inter)
    return(rawp)
  }
  
  # parallel simulation
  if (!is.null(file_cache) && file.exists(as.character(file_cache))){
    df_null <- read_rds(file_cache)
  } else {
    # set parallel processing
    Ns_iter <- rep(N_subj, times=iter)
    ls_tibble <- pbapply::pblapply(Ns_iter, sim_hber_single, cl=N_core)
    
    df_null <- reduce(ls_tibble, rbind) %>% 
      dplyr::mutate(iter=1:n())
    
    if (!is.null(file_cache)){
      write_rds(df_null, file=file_cache)
    }
  }
  
  return(df_null)
}


### Apply alpha to df
sim_train_cber <- function (df_sim, alpha=.05, disp=TRUE) {
  # function to calculate the Type I error rate. 
  
  # df_sim: dataframe obtained from sim_train_null
  # alpha: alpha level used to calculate the Type I error rate.
  # disp: if true, the output will be a wide format. If false, it will be long.
  
  df_sig <- df_sim %>% 
    mutate(alpha = alpha,
           sig_d_simple = p_d_simple <= alpha, # apply alpha
           sig_d_inter = p_d_inter <= alpha, # apply alpha
           sig_rt_simple = p_rt_simple <= alpha, # apply alpha
           sig_rt_inter = p_rt_inter <= alpha, # apply alpha
           sig_d_both = sig_d_simple * sig_d_inter, # whether both are sig for d
           sig_rt_both = sig_rt_simple * sig_rt_inter, # whether both are sig for rt
           sig_d_or_rt = sig_d_both | sig_rt_both) %>% # EOI
    select(N_subj, starts_with("sig_")) %>% 
    pivot_longer(starts_with("sig_"), names_to = "effects", values_to = "isSig") %>% 
    group_by(N_subj, effects) %>% 
    summarize(N_iter = n(),
              N_sig = sum(isSig),
              Type_I = mean(isSig),
              .groups="keep")
  
  if (disp){
    df_sig <- df_sig %>% 
      select(-c(N_sig)) %>% 
      pivot_wider(c(N_subj, N_iter), names_from = "effects", values_from = "Type_I")
  }
  
  return(df_sig)
}
