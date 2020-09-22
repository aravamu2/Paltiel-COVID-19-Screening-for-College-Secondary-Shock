library(tidyverse)
library(here)

### function to run model based on input parameters -------------------------------
run_model <- function(exogenous.shocks,
                      frequency.exogenous.shocks.per.day,
                      frequency.screening,
                      time.to.recovery,
                      per.asymptotics.advancing.to.symptoms,
                      r0,
                      symptom.case.fatality.ratio,
                      days.incubation,
                      time.to.return.fps.from.isolation,
                      new.frequency.screening,
                      new.r0,
                      initial.susceptible,
                      initial.infected,
                      new.infections.per.shock,
                      secondary.shock,
                      new.day,
                      new.infections,
                      test.specificity,
                      test.sensitivity,
                      test.cost = NULL,
                      confirmatory.test.cost = NULL) {
  
  ## calculated parameters 
  num.exogenous.shocks <- case_when(
    exogenous.shocks == "Yes" ~ 1,
    exogenous.shocks == "No" ~ 0
  )
  cycles.per.day <- 3
  frequency.exogenous.shocks <- cycles.per.day*frequency.exogenous.shocks.per.day
  cycles.per.test <- case_when(
    frequency.screening == "Daily" ~ 1*cycles.per.day,
    frequency.screening == "Every 2 Days" ~ 2*cycles.per.day,
    frequency.screening == "Every 3 Days" ~ 3*cycles.per.day,
    frequency.screening == "Weekly" ~ 7*cycles.per.day,
    frequency.screening == "Every 2 Weeks" ~ 14*cycles.per.day,
    frequency.screening == "Every 3 Weeks" ~ 21*cycles.per.day,
    frequency.screening == "Every 4 Weeks" ~ 28*cycles.per.day,
    frequency.screening == "Symptoms Only" ~ 99999999999
  )
  rho <- 1/(time.to.recovery*cycles.per.day)
  sigma <- rho*(per.asymptotics.advancing.to.symptoms/(1-per.asymptotics.advancing.to.symptoms))
  beta <- r0*(rho+sigma)
  delta <- (symptom.case.fatality.ratio/(1-symptom.case.fatality.ratio))*rho
  theta <- 1/(days.incubation*cycles.per.day)
  mu <- 1/(cycles.per.day*time.to.return.fps.from.isolation)
  
  new.cycles.per.test <- case_when(
    new.frequency.screening == "Daily" ~ 1*cycles.per.day,
    new.frequency.screening == "Every 2 Days" ~ 2*cycles.per.day,
    new.frequency.screening == "Every 3 Days" ~ 3*cycles.per.day,
    new.frequency.screening == "Weekly" ~ 7*cycles.per.day,
    new.frequency.screening == "Every 2 Weeks" ~ 14*cycles.per.day,
    new.frequency.screening == "Every 3 Weeks" ~ 21*cycles.per.day,
    new.frequency.screening == "Every 4 Weeks" ~ 28*cycles.per.day,
    new.frequency.screening == "Symptoms Only" ~ 99999999999
  )
  new.beta <- new.r0*(rho+sigma)
  
  
  ## calculate compartment matrix 
  n.cycle <- 240
  
  mat <- matrix(c(0,initial.susceptible,0,0,initial.infected,0,0,0,0), nrow = 1)
  mat <- rbind(
    mat,
    c(1,
      max(0,mat[1,2]*(1-beta*(mat[1,5]/(mat[1,2]+mat[1,5]+mat[1,4])))+mat[1,3]*mu),
      max(0,mat[1,3]*(1-mu)),
      max(0,mat[1,4]*(1-theta)+ beta*(mat[1,2]*mat[1,5]/(mat[1,2]+mat[1,5]+mat[1,4]))),
      max(0,mat[1,5]*(1-sigma-rho)+mat[1,4]*theta),
      max(0,mat[1,6]*(1-delta-rho)+(mat[1,5]+mat[1,7])*sigma),
      0,
      max(0,mat[1,8]+(mat[1,5]+mat[1,6]+mat[1,7])*rho),
      max(0,delta*mat[1,6]+mat[1,9]))
  )
  
  superspreader.event <- 0
  superspreader.event <- c(superspreader.event, 
                           (1:n.cycle %% frequency.exogenous.shocks == 0)*num.exogenous.shocks)
  superspreader.infections <- superspreader.event * new.infections.per.shock
  if(secondary.shock == "Yes") {
    superspreader.event[new.day*3] <- 1
    superspreader.infections[new.day*3+1] <- superspreader.infections[new.day*3+1] + new.infections
  }
  
  for(i in 2:n.cycle) {
    if((secondary.shock == "No") | (new.day > (i/3))) {
      mat <- 
        rbind(mat,
              c(i,
                max(0,mat[i,2]*(1-beta*(mat[i,5]/(mat[i,2]+mat[i,5]+mat[i,4])))+mat[i,3]*mu-mat[i-1,2]*(1-test.specificity)/cycles.per.test-superspreader.event[i+1]*superspreader.infections[i+1]),
                max(0,mat[i,3]*(1-mu)+mat[i-1,2]*(1-test.specificity)/cycles.per.test),
                max(0,mat[i,4]*(1-theta)+beta*(mat[i,2]*mat[i,5]/(mat[i,2]+mat[i,5]+mat[i,4]))+superspreader.event[i+1]*superspreader.infections[i+1]),
                max(0,mat[i,5]*(1-sigma-rho)+mat[i,4]*theta-mat[i-1,5]*test.sensitivity/cycles.per.test),
                max(0,mat[i,6]*(1-delta-rho)+(mat[i,5]+mat[i,7])*sigma),
                max(0,mat[i,7]*(1-sigma-rho)+mat[i-1,5]*test.sensitivity/cycles.per.test),
                max(0,mat[i,8]+(mat[i,5]+mat[i,6]+mat[i,7])*rho),
                max(0,delta*mat[i,6]+mat[i,9]))
        )
    } else if((secondary.shock == "Yes") & (new.day <= (i/3))) {
      mat <- 
        rbind(
          mat,
          c(i,
            max(0,mat[i,2]*(1-new.beta*(mat[i,5]/(mat[i,2]+mat[i,5]+mat[i,4])))+mat[i,3]*mu-mat[i-1,2]*(1-test.specificity)/new.cycles.per.test-superspreader.event[i+1]*superspreader.infections[i+1]),
            max(0,mat[i,3]*(1-mu)+mat[i-1,2]*(1-test.specificity)/new.cycles.per.test),
            max(0,mat[i,4]*(1-theta)+new.beta*(mat[i,2]*mat[i,5]/(mat[i,2]+mat[i,5]+mat[i,4]))+superspreader.event[i+1]*superspreader.infections[i+1]),
            max(0,mat[i,5]*(1-sigma-rho)+mat[i,4]*theta-mat[i-1,5]*test.sensitivity/new.cycles.per.test),
            max(0,mat[i,6]*(1-delta-rho)+(mat[i,5]+mat[i,7])*sigma),
            max(0,mat[i,7]*(1-sigma-rho)+mat[i-1,5]*test.sensitivity/new.cycles.per.test),
            max(0,mat[i,8]+(mat[i,5]+mat[i,6]+mat[i,7])*rho),
            max(0,delta*mat[i,6]+mat[i,9]))
        )
    }
  }
  
  mat <- cbind(mat, superspreader.event, superspreader.infections)
  
  names.df <- c("Cycle","Susceptible","FP","Exposed","Asympt","Symptoms","TP","Recovered","Dead","Superspreader Event", "Superspreader Infections")
  df <- 
    mat %>% 
    as_tibble() %>% 
    rename_all(~names.df) %>% 
    mutate(`Persons Tested` = ifelse(secondary.shock == "No" | new.day > ((Cycle+1)/3),(lag(Susceptible,1,NA)+lag(Exposed,1,NA)+lag(Asympt,1,NA))/cycles.per.test,(lag(Susceptible,1,NA)+lag(Exposed,1,NA)+lag(Asympt,1,NA))/new.cycles.per.test),
           `Total TPs` = ifelse(secondary.shock == "No" | new.day > ((Cycle+2)/3),lag(Asympt,2,NA)*test.sensitivity/cycles.per.test,lag(Asympt,2,NA)*test.sensitivity/new.cycles.per.test),
           `Total FPs` = ifelse(secondary.shock == "No" | new.day > ((Cycle+2)/3),lag(Susceptible,2,NA)*(1-test.specificity)/cycles.per.test,lag(Susceptible,2,NA)*(1-test.specificity)/new.cycles.per.test),
           `Total TNs` = ifelse(secondary.shock == "No" | new.day > ((Cycle+2)/3),lag(Susceptible,2,NA)*test.specificity/cycles.per.test,lag(Susceptible,2,NA)*test.specificity/new.cycles.per.test),
           `Total FNs` = ifelse(secondary.shock == "No" | new.day > ((Cycle+2)/3),lag(Exposed,2,NA)+lag(Asympt,2,NA)*(1-test.sensitivity)/cycles.per.test,lag(Exposed,2,NA)+lag(Asympt,2,NA)*(1-test.sensitivity)/new.cycles.per.test)) %>% 
    mutate(Day = Cycle/cycles.per.day,
           `True Positive` = TP,
           Symptoms = Symptoms,
           `False Positive` = FP,
           Total = TP+Symptoms+FP) %>% 
    mutate(`New Infections` = ifelse(secondary.shock == "No" | new.day > ((Cycle+1)/3),lag(Asympt,1,NA)*beta*lag(Susceptible,1,NA)/(lag(Susceptible,1,NA)+lag(Exposed,1,NA)+lag(Asympt,1,NA)),lag(Asympt,1,NA)*new.beta*lag(Susceptible,1,NA)/(lag(Susceptible,1,NA)+lag(Exposed,1,NA)+lag(Asympt,1,NA))),
           `New Infections` = ifelse(Cycle>1,
                                     `New Infections`+pmin(`Superspreader Infections`,lag(Susceptible,1,NA)),
                                     `New Infections`),
           `New Infections` = ifelse(is.na(`New Infections`),0,`New Infections`),
           `Cumulative Infections` = cumsum(`New Infections`),
           `%Cumulative Infections` = `Cumulative Infections`/initial.susceptible)
  
}

### summary statistics function ------------------------------------------------

summary_stats <- function(df) {
  # df is output from `run_model()`
  df %>% 
    slice(2:n()) %>% 
    summarize(`Peak Isolation Occupancy` = max(`Total`, na.rm = TRUE),
              `Time of Peak Isolation` = `Day`[which.max(`Total`)],
              `Total Persons Tested in 80 days` = sum(`Persons Tested`, na.rm = TRUE),
              `Total Confirmatory Tests Performed` = sum(`Total TPs`, na.rm = TRUE) + sum(`Total FPs`, na.rm = TRUE),
              `Average Isolation Unit Census` = mean(`Total`, na.rm = TRUE),
              `Average %TP in Isolation` = 1-(mean(`False Positive`, na.rm = TRUE)/mean(`Total`, na.rm = TRUE)),
              `Total testing cost` = `Total Persons Tested in 80 days`*test.cost+`Total Confirmatory Tests Performed`*confirmatory.test.cost,
              `Total Infections` = last(`Cumulative Infections`))
  
}

### input parameters for sweep -------------------------------------------------

param_sweep <- expand_grid(
  initial.susceptible = 6800,
  initial.infected = c(outer(c(1,5),10^(1:2),"*"),1000),
  
  r0 = seq(0.5,1.5,0.25),
  exogenous.shocks = "Yes",
  frequency.exogenous.shocks.per.day = 7,
  new.infections.per.shock = c(outer(c(1,5),10^(0:1),"*"),100),
  
  days.incubation = 3,
  time.to.recovery = 14,
  per.asymptotics.advancing.to.symptoms = 0.3,
  symptom.case.fatality.ratio = 0.0005,
  
  frequency.screening = c("Symptoms Only",
                          "Every 4 Weeks",
                          "Every 3 Weeks",
                          "Every 2 Weeks",
                          "Weekly",
                          "Every 3 Days",
                          "Every 2 Days",
                          "Daily"),
  test.sensitivity = 0.8,
  test.specificity = 0.98,
  test.cost = 25,
  time.to.return.fps.from.isolation = c(1,2,14),
  confirmatory.test.cost = 100,
  
  secondary.shock = "No", # no secondary shock for now, ignore below variables
  new.day = 40,
  new.infections = 50,
  new.frequency.screening = "Weekly",
  new.r0 = 2.5
)

param_sweep <-
  param_sweep %>% 
  rowid_to_column() %>% 
  group_by(rowid)  %>% 
  nest(params = -rowid)

### run everything  ------------------------------------------------------------

res <-
  param_sweep %>% 
  mutate(df = map(params, ~do.call(run_model, args = .x))) %>% 
  mutate(sum_stat = map(df, summary_stats)) %>% 
  select(-df) %>% 
  unnest_wider(params) %>% 
  unnest_wider(sum_stat) 

res

date <- Sys.Date() %>% str_remove_all("-")
write_csv(res, path = here(paste0("output_parameter_sweep_", date, ".csv")))
  

