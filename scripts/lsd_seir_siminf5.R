library(SimInf)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

##################################################################
## MODEL INITIALIZATION

# define compartments
compartments <- c("Sc", "Ec1", "Ec2asymp", "Ec2symp", "Icsymp", "Icasymp", 
                  "Icpredeath", "Icprerecovery", "Icum", "Isympcum", "Rc", "Dc", "Culled",
                  "Sv", "Iv", "Ivcum")

# define transitions between compartments
transitions <- c(#cattle
  "Sc -> Nc>0 ? k*p_c*Sc*Iv/Nc : 0 -> Ec1", 
  "Ec1 -> a*ps*Ec1 -> Ec2symp+ Icum + Isympcum",  "Ec1 -> a*(1-ps)*Ec1 -> Ec2asymp + Icum",
  "Ec2asymp -> gamma*Ec2asymp -> Icasymp ", 
  "Icasymp -> tau_CA*Icasymp -> Rc",
  "Ec2symp -> gamma*Ec2symp -> Icsymp ", 
  "Icsymp -> Icsymp*b*pd -> Icpredeath" , "Icsymp -> Icsymp*b*(1-pd) -> Icprerecovery",
  "Icpredeath -> alpha*Icpredeath -> Dc",
  "Icprerecovery -> tau_CS*Icprerecovery -> Rc", 
  # vectors
  "Sv -> Nc>0 ? (k*Sv/Nc)*(p_vs*(Icsymp+Icprerecovery+Icpredeath)+p_va*Icasymp) : 0 -> Iv + Ivcum",
  "Iv -> tau_v*Iv -> Sv", 
  "@ -> Sv*mu -> Sv", "Sv -> mu*Sv -> @", "Iv -> mu*Iv -> @" ,
  # variables
  "Nc <- Sc+Ec1+Ec2asymp+Ec2symp+Icsymp+Icasymp+Icpredeath+Icprerecovery+Rc"
)


# initial conditions
n <- 100 #number of stochastic repetitions
u0 <- data.frame(Sc = rep(99, n), Ec1 = rep(1, n), Ec2asymp = rep(0, n), Ec2symp = rep(0, n),
                 Icsymp = rep(0, n), Icasymp = rep(0, n), Icpredeath = rep(0, n), Icprerecovery = rep(0, n),
                 Icum = rep(0, n), Isympcum = rep(0,n),
                 Rc = rep(0, n), Dc = rep(0,n),
                 Culled=rep(0,n),
                 Sv=rep(5000,n), Iv=rep(0,n), Ivcum=rep(0,n))


#######################################################################
# INITIALIZE SCHEDULED EVENTS 

n_events<-7
E <- matrix(rep(0,ncol(u0)*n_events), nrow = ncol(u0), ncol = n_events, 
            dimnames = list(names(u0),as.character(1:n_events))) 
N <- matrix(rep(0,ncol(u0)*n_events), nrow = ncol(u0), ncol = n_events, 
            dimnames = list(names(u0),as.character(1:n_events))) 

E["Ec2symp","1"] <- 1 #sample from Ec2symp compartment
N["Ec2symp","1"] <- 9 #move them 7 compartments forward, to culled

E["Ec2asymp","2"] <- 1 #sample from Ec2asymp compartment
N["Ec2asymp","2"] <- 10 #move them 8 compartments forward, to culled

E["Icsymp","3"] <- 1 #sample from Icsymp compartment
N["Icsymp","3"] <- 8 #move them 6 compartments forward, to culled


E["Icpredeath","4"] <- 1 
N["Icpredeath","4"] <- 6 

E["Icprerecovery","5"] <- 1 
N["Icprerecovery","5"] <- 5 

E["Icasymp","6"] <- 1 
N["Icasymp","6"] <- 7 

E["Sv","7"] <- 1 
E["Iv","7"] <- 1


#######################################################################
# TEST DIFFERENT STRATEGIES

# tested strategies
days_detectable <- seq(1,10,1)
se_asymp <- seq(0,1,0.1)
effectiveness_insecticide <- c(0, 0.5, 0.8)

test_combinations <- expand.grid(days_detectable, se_asymp, effectiveness_insecticide)
names(test_combinations) <- c("days_detectable","se_asymp", "effectiveness_insecticide")

iip <- 10
r <- 0.5 #reduction se for latent stage 2 animals
start_interventions <- 7

#initialize result dfs
n_metrics <- 5

size_results <- nrow(test_combinations) * n_metrics
results <- data.frame(
  se = numeric(size_results),
  days_detectable = numeric(size_results),
  effectiveness_insecticide = numeric(size_results),
  metric = numeric(size_results),
  med = numeric(size_results),
  avg = numeric(size_results),
  p5 = numeric(size_results),
  p95 = numeric(size_results)
)

size_df_log <- n*nrow(test_combinations)
df_log <- data.frame(
  node = numeric(size_df_log),
  tot_cases_cattle = numeric(size_df_log),
  tot_cases_vector = numeric(size_df_log),
  tot_culled = numeric(size_df_log),
  duration_epi = numeric(size_df_log),
  duration_epi_tot = numeric(size_df_log),
  se = numeric(size_df_log),
  days_detectable = numeric(size_df_log),
  effectiveness_insecticide = numeric(size_df_log)
)

row_id_res <- 1
row_id_log <- 1

for (i in 1:nrow(test_combinations)) {
  
  cur_days_detectable <- test_combinations[i, "days_detectable"]
  cur_se <- test_combinations[i, "se_asymp"]
  cur_eff_insecticide <- test_combinations[i, "effectiveness_insecticide"]
  
  if(cur_days_detectable == iip) {
    a_dur=1e-5
  } else {
    a_dur=iip-cur_days_detectable
  }
  
  #model parameters
  params <- c(k=0.5, p_c=0.07,  
              a=1/a_dur, ps=0.5,
              gamma=1/cur_days_detectable, 
              b=1, pd = 0.12,
              tau_CA=1/15,tau_CS=1/40, alpha=1/20,
              #p_vs=0.22, p_va=0.006,
              p_vs=0.46, p_va=0.2,
              nu=1/21,
              tau_v=1/5.5, mu=1/21)
  
  #run model init to get initial conditions for simu
  model_init <- mparse(transitions = transitions, compartments = compartments,
                       #parameter values
                       gdata = params, 
                       u0 = u0, tspan = 1:200)
  
  # run simulations 
  result_init <- run(model = model_init)
  traj_init<-trajectory(result_init)
  
  traj_init_long <- traj_init %>% rename("incidence"="Icum", "incidence_symp"="Isympcum", "incidence_v"="Ivcum") %>% 
    select(contains("Iv"), contains("Ic"),Ec1,Ec2asymp,Ec2symp, Dc, node, time, incidence, incidence_symp, incidence_v) %>% 
    group_by(node) %>%
    mutate(incidence = diff(c(0, incidence)), 
           incidence_symp = diff(c(0, incidence_symp)), incidence_v = diff(c(0, incidence_v)),
           infected=Ec2asymp+Ec2symp+Icsymp+Icasymp+Icpredeath+Icprerecovery,
           infected_sympto=Ec2symp+Icsymp+Icpredeath+Icprerecovery) %>%
    select(infected_sympto, infected, Dc, time, node) %>%
    pivot_longer(cols = -c(node, time), names_to = "compartment",values_to = "count")
  
  # plot and save trajectories for inspection
  # traj_init_stats <- traj_init_long %>%
  #   group_by(time, compartment) %>%
  #   summarise(med = median(count),
  #             q05  = quantile(count, 0.05), q95  = quantile(count, 0.95),
  #             .groups = "drop")
  # t_plot<-ggplot(traj_init_stats , aes(x = time, color = compartment, fill = compartment)) +
  #   geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2, colour = NA) +
  #   geom_line(aes(y = med), linewidth = 1) +
  #   theme_minimal() +
  #   labs(x = "Time (days)",y = "Count",title = "Example of trajectories of the model with no interventions") +
  #   theme_bw() +
  #   scale_color_manual(name="",
  #                      labels = c("Cumulative deaths", "Infected", "Symptomatic infected"),
  #                      values = c("#5A9CB5", "#FAAC68", "#FA6868"))+
  #   scale_fill_manual(name="",
  #                     labels = c("Cumulative deaths", "Infected", "Symptomatic infected"),
  #                     values = c("#5A9CB5", "#FAAC68", "#FA6868"))+
  #   theme(legend.position = c(0.95, 0.95),
  #     legend.justification = c(1, 1),
  #     legend.background = element_rect(fill = "white", color = "black"),
  #     legend.title = element_blank())

  # extract time with first symptomatic
  time_detection <- traj_init %>% filter(Icpredeath>0 | Icprerecovery>0 | Icsymp>0)  %>% 
    group_by(node) %>% summarise(time_detection=min(time))
  
  # inspection of the time to detection
  # d_plot<-ggplot(time_detection) + geom_histogram(aes(x=time_detection), fill="#9A3F3F", col="black") + theme_bw() +
  #   labs(x="Detection time (first symptomatic animal)", y="Count") #title="Distribution of the timing of detection over 100 simulations"
  # ggsave(paste0("figures/final_param_set/time_detection/time_detection_",i,".png"), bg="white")
  # ggarrange(t_plot,d_plot, ncol = 1, nrow = 2, labels=c('A','B'))
  # ggsave("figures/final_param_set/fig2_2panels.png", bg="white")
  
  # create new initial conditions to start new simulations
  new_init_values <- traj_init %>%
    semi_join(time_detection, by = c("node" = "node", "time" = "time_detection")) %>%
    arrange(node) %>%
    select(-c(node,time))
  
  # define control strategies as events
  test_cull_Ecsymp_cur <- data.frame(event = "intTrans", #internal transfer
                                     time=rep(seq(start_interventions,500, by=2), each=nrow(new_init_values)), #every second day
                                     node= 1:nrow(new_init_values), #this indicates that it should be applied for all stochastic repetitions
                                     dest = 0, n = 0,
                                     proportion = r, #probability of detection
                                     select = 1, shift = 1) #columns in E and N to be used
  
  test_cull_Ecasymp_cur <- data.frame(event = "intTrans", 
                                      time=rep(seq(start_interventions,500, by=2), each=nrow(new_init_values)), 
                                      node= 1:nrow(new_init_values), 
                                      dest = 0, n = 0,
                                      proportion = cur_se, 
                                      select = 2, shift = 2)
  
  test_cull_Isymp_cur <- data.frame(event = "intTrans", 
                                    time=rep(seq(start_interventions,500, by=2), each=nrow(new_init_values)), 
                                    node= 1:nrow(new_init_values), 
                                    dest = 0, n = 0,
                                    proportion = 1, 
                                    select = 3, shift = 3)
  
  test_cull_predeath_cur <- data.frame(event = "intTrans",
                                       time=rep(seq(start_interventions,500, by=2), each=nrow(new_init_values)),
                                       node= 1:nrow(new_init_values),
                                       dest = 0, n = 0,
                                       proportion = 1,
                                       select = 4, shift = 4)
  
  test_cull_prerecov_cur <- data.frame(event = "intTrans",
                                       time=rep(seq(start_interventions,500, by=2), each=nrow(new_init_values)),
                                       node= 1:nrow(new_init_values),
                                       dest = 0, n = 0,
                                       proportion = 1,
                                       select = 5, shift = 5)
  
  test_cull_Iasymp_cur <- data.frame(event = "intTrans", 
                                     time=rep(seq(start_interventions,500, by=2), each=nrow(new_init_values)), 
                                     node= 1:nrow(new_init_values), 
                                     dest = 0, n = 0,
                                     proportion = cur_se, 
                                     select = 6, shift = 6)
  
  insecticide <- data.frame(event = "exit", 
                            time=rep(start_interventions, each=nrow(new_init_values)), 
                            node= 1:nrow(new_init_values), 
                            dest = 0, n = 0,
                            proportion = cur_eff_insecticide, 
                            select = 7, shift = 0)
  
  all_events_cur<-rbind(test_cull_Ecsymp_cur, test_cull_Isymp_cur, test_cull_predeath_cur,
                        test_cull_prerecov_cur, test_cull_Ecasymp_cur, test_cull_Iasymp_cur, insecticide)
  
  # run model with events starting from initial values defined previously
  model_events <- mparse(transitions = transitions, compartments = compartments,
                         #parameter values
                         gdata = params, 
                         u0 = new_init_values, tspan = 1:500,
                         events=all_events_cur,
                         E=E, N=N)
  
  result_final <- run(model = model_events)
  
  traj_final <- trajectory(result_final)
  
  traj_final_long <- traj_final %>% rename("incidence"="Icum", "incidence_symp"="Isympcum", "incidence_v"="Ivcum") %>% 
    select(contains("Iv"), contains("Ic"), Culled, node, time, incidence, incidence_symp, incidence_v) %>% 
    group_by(node) %>%
    mutate(Ic_tot = rowSums(across(contains("Ic"))), incidence = diff(c(0, incidence)), 
           incidence_symp = diff(c(0, incidence_symp)), incidence_v = diff(c(0, incidence_v)),
           Ic_sympto=Icsymp+Icpredeath+Icprerecovery) %>%
    pivot_longer(cols = -c(node, time), names_to = "compartment",values_to = "count")
  
  # plot and save trajectories for inspection
  # p_t<-ggplot(traj_final_long %>% filter(compartment %in% c("Iv","Ic_tot","Ic_sympto","Culled"), node %in% 1:50),
  #             aes(x=time, y=count, group=node, color=compartment)) +
  #   geom_line() + ggtitle("Trajectory of the model with test and remove strategy") +
  #   facet_wrap(~compartment, scale="free_y",
  #              labeller = as_labeller(c(Ic_sympto = "Symptomatic animal",
  #                                       Ic_tot = "Infected animals", Culled="Culled animal",
  #                                       Iv = "Infected vector")))  +
  #   theme_bw() + theme(legend.position = "none") +
  #   scale_color_manual(values = c("#5A9CB5", "#FACE68", "#FAAC68", "#FA6868")) + xlim(0,150)
  # ggsave(paste0("figures/final_param_set/trajectories_se_daysd_500insect_2/ex_trajectory_curdaysd",cur_days_detectable,"_se",cur_se,"_insect",cur_eff_insecticide,".png"), bg="white")
  
  duration_epi <- traj_final_long %>% filter(compartment=="Ic_tot", count>0) %>% 
    group_by(node) %>% summarise(duration_epi=max(time)) %>% 
    left_join(time_detection) %>%
    mutate(duration_epi_tot=time_detection+duration_epi) %>% select(-time_detection) 
  # pivot_longer(-node, names_to="type", values_to = "time") #for plotting
  
  # ggplot(duration_epi2 %>% filter(type!="duration_epi_tot")) + geom_line(aes(x=type, y=time, group=node), col="black") + theme_bw() +
  #    labs(x="Detection time (first symptomatic animal)", y="Count") 
  
  metrics_all_nodes <- traj_final_long %>% 
    pivot_wider(names_from = compartment, values_from = count) %>%
    group_by(node) %>% 
    summarise(tot_cases_cattle = sum(incidence), tot_cases_vector=sum(incidence_v), tot_culled=max(Culled)) %>% 
    ungroup() %>% left_join(duration_epi) %>%
    mutate(se=rep(cur_se, nrow(.)), days_detectable=rep(cur_days_detectable, nrow(.)), 
           effectiveness_insecticide=rep(cur_eff_insecticide, nrow(.)))
  df_log[row_id_log:(row_id_log + nrow(metrics_all_nodes) - 1), ] <- metrics_all_nodes
  
  metrics_summary <- metrics_all_nodes %>% select(-c(se,days_detectable, effectiveness_insecticide)) %>%
    summarise(med_duration=median(duration_epi, na.rm = TRUE), 
              avg_duration=mean(duration_epi, na.rm = TRUE),
              p5_duration=quantile(duration_epi, 0.05, na.rm = TRUE), 
              p95_duration=quantile(duration_epi, 0.95, na.rm = TRUE),
              med_durationtot=median(duration_epi_tot, na.rm = TRUE), 
              avg_durationtot=mean(duration_epi_tot, na.rm = TRUE),
              p5_durationtot=quantile(duration_epi_tot, 0.05, na.rm = TRUE), 
              p95_durationtot=quantile(duration_epi_tot, 0.95, na.rm = TRUE),
              med_Ic=median(tot_cases_cattle), 
              avg_Ic=mean(tot_cases_cattle), 
              p5_Ic=quantile(tot_cases_cattle, 0.05, na.rm = TRUE), 
              p95_Ic=quantile(tot_cases_cattle, 0.95, na.rm = TRUE),
              med_culled=median(tot_culled), 
              avg_culled=mean(tot_culled), 
              p5_culled=quantile(tot_culled, 0.05, na.rm = TRUE), 
              p95_culled=quantile(tot_culled, 0.95, na.rm = TRUE),
              med_Iv=median(tot_cases_vector),
              avg_Iv=mean(tot_cases_vector),
              p5_Iv=quantile(tot_cases_vector, 0.05, na.rm = TRUE), 
              p95_Iv=quantile(tot_cases_vector, 0.95, na.rm = TRUE)) %>%
    mutate(se=rep(cur_se, nrow(.)), days_detectable=rep(cur_days_detectable, nrow(.)), 
           effectiveness_insecticide=rep(cur_eff_insecticide, nrow(.))) %>%
    pivot_longer(cols=-c(se, days_detectable,effectiveness_insecticide), 
                 names_to = c("stat","metric"), names_sep = "_", values_to = "value") %>%
    pivot_wider(names_from = stat,values_from = value)
  
  results[row_id_res:(row_id_res + n_metrics - 1), ] <- metrics_summary
  
  row_id_res <- row_id_res + n_metrics
  row_id_log <- row_id_log + n
}


# save(results, file = "Results/result_summary_metric_100rep.RData")
# save(df_log, file = "Results/result_log_100rep.RData")

#######################################################################
# PLOT RESULTS

load("Results/result_summary_metric_100rep.RData")
load("Results/result_log_100rep.RData")

global_limits_culled <- range((results %>% filter(metric=="culled"))$avg, na.rm = TRUE)
global_limits_dur <- range((results %>% filter(metric=="duration"))$avg, na.rm = TRUE)
global_limits_ic <- range((results %>% filter(metric=="Ic"))$avg, na.rm = TRUE)
global_limits_iv <- range((results %>% filter(metric=="Iv"))$avg, na.rm = TRUE)

for (i in effectiveness_insecticide) { 
  
  
  p1<-ggplot(results %>% filter(metric=="culled", effectiveness_insecticide==i), aes(x = se, y = days_detectable, fill = avg)) +
    geom_tile() +
    scale_fill_viridis_c(limits=global_limits_culled) +
    theme_bw() +
    labs(
      x = "Sensitivity of diagnostic test for asymptomatics",
      y = "Number of days when the animal\nis detectable before being infectious",
      fill = "Value",
      title="Number of culled animals"
    )
  
  p2<-ggplot(results %>% filter(metric=="duration", effectiveness_insecticide==i), aes(x = se, y = days_detectable, fill = avg)) +
    geom_tile() +
    scale_fill_viridis_c(limits=global_limits_dur) +
    #geom_text(aes(label = round(avg, 1)), size = 3) +
    theme_bw() +
    labs(
      x = "Sensitivity of diagnostic test for asymptomatics",
      y = "Number of days when the animal\nis detectable before being infectious",
      fill = "Value",
      title="Time to virus extinction"
    )
  
  p3<-ggplot(results %>% filter(metric=="Iv", effectiveness_insecticide==i), aes(x = se, y = days_detectable, fill = avg)) +
    geom_tile() +
    scale_fill_viridis_c(limits=global_limits_iv) +
    theme_bw() +
    labs(
      x = "Sensitivity of diagnostic test for asymptomatics",
      y = "Number of days when the animal\nis detectable before being infectious",
      fill = "Value",
      title="Total number of\ninfected vectors"
    )
  
  p4<-ggplot(results %>% filter(metric=="Ic", effectiveness_insecticide==i), aes(x = se, y = days_detectable, fill = avg)) +
    geom_tile() +
    scale_fill_viridis_c(limits=global_limits_ic) +
    #geom_text(aes(label = round(avg, 1)), size = 3) +
    theme_bw() +
    labs(
      x = "Sensitivity of diagnostic test for asymptomatics",
      y = "Number of days when the animal\nis detectable before being infectious",
      fill = "Value",
      title="Total number of\ninfected cows"
    )
  
  ggarrange(p1,p2,p3,p4, nrow=2, ncol=2)
  ggarrange(
    p1 + theme(axis.title= element_blank()),
    p2 + theme(axis.title = element_blank()),
    p3 + theme(axis.title = element_blank()),
    p4 + theme(axis.title = element_blank()),
    nrow = 2, ncol = 2,
    align = "hv", labels=c('A','B','C','D')) %>%
    annotate_figure(
      left = text_grob("Number of days when the animal\nis detectable before being infectious",  rot = 90),
      bottom = text_grob("Sensitivity of diagnostic test for asymptomatics")
    )
  ggsave(paste0("figures/final_param_set/heatmap5_insecteff",i,"_avg.png"), bg="white")
}


# summarise results

results %>% group_by(effectiveness_insecticide, metric) %>% summarise(min=min(avg), med=median(avg), max=max(avg)) %>% arrange(metric)

a1 <- ggplot(results %>% filter(metric=="duration"), aes(x = se, y = days_detectable, fill = avg)) +
  geom_tile() +
  scale_fill_viridis_c() +
  # geom_text(aes(label = round(avg, 1)), size = 3) +
  facet_wrap(~effectiveness_insecticide) +
  theme_bw() +
  labs(
    x = "Sensitivity of diagnostic test for asymptomatics",
    y = "Number of days when the animal\nis detectable before being infectious",
    fill = "Value",
    title="Time to virus extinction"
  )

a2 <- ggplot(results %>% filter(metric=="Ic"), aes(x = se, y = days_detectable, fill = avg)) +
  geom_tile() +
  scale_fill_viridis_c() +
  # geom_text(aes(label = round(avg, 1)), size = 3) +
  facet_wrap(~effectiveness_insecticide) +
  theme_bw() +
  labs(
    x = "Sensitivity of diagnostic test for asymptomatics",
    y = "Number of days when the animal\nis detectable before being infectious",
    fill = "Value",
    title="Total number of infected cows"
  )

ggarrange(
  a1 + theme(axis.title= element_blank()),
  a2 + theme(axis.title = element_blank()),
  nrow = 2, ncol = 1,
  align = "hv", labels=c('A','B')) %>%
  annotate_figure(
    left = text_grob("Number of days when the animal\nis detectable before being infectious",  rot = 90),
    bottom = text_grob("Sensitivity of diagnostic test for asymptomatics")
  )
ggsave("figures/final_param_set/heatmap5_3diffinsecteff_avg.png", bg="white")



p1<-ggplot(results %>% filter(metric=="culled", effectiveness_insecticide==i), aes(x = se, y = days_detectable, fill = p95)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(
    x = "Sensitivity of diagnostic test for asymptomatics",
    y = "Number of days when the animal\nis detectable before being infectious",
    fill = "Value",
    title="Number of culled animals (95th percentile)"
  )

p2<-ggplot(results %>% filter(metric=="duration", effectiveness_insecticide==i), aes(x = se, y = days_detectable, fill = p95)) +
  geom_tile() +
  scale_fill_viridis_c() +
  #geom_text(aes(label = round(avg, 1)), size = 3) +
  theme_bw() +
  labs(
    x = "Sensitivity of diagnostic test for asymptomatics",
    y = "Number of days when the animal\nis detectable before being infectious",
    fill = "Value",
    title="Time to virus extinction (95th percentile)"
  )

p3<-ggplot(results %>% filter(metric=="Iv", effectiveness_insecticide==i), aes(x = se, y = days_detectable, fill = p95)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(
    x = "Sensitivity of diagnostic test for asymptomatics",
    y = "Number of days when the animal\nis detectable before being infectious",
    fill = "Value",
    title="Total number of\ninfected vectors (95th percentile)"
  )

p4<-ggplot(results %>% filter(metric=="Ic", effectiveness_insecticide==i), aes(x = se, y = days_detectable, fill = p95)) +
  geom_tile() +
  scale_fill_viridis_c() +
  #geom_text(aes(label = round(avg, 1)), size = 3) +
  theme_bw() +
  labs(
    x = "Sensitivity of diagnostic test for asymptomatics",
    y = "Number of days when the animal\nis detectable before being infectious",
    fill = "Value",
    title="Total number of\ninfected cows (95th percentile)"
  )

ggarrange(p1,p2,p3,p4, nrow=2, ncol=2)
ggarrange(
  p1 + theme(axis.title= element_blank()),
  p2 + theme(axis.title = element_blank()),
  p3 + theme(axis.title = element_blank()),
  p4 + theme(axis.title = element_blank()),
  nrow = 2, ncol = 2,
  align = "hv", labels=c('A','B','C','D')) %>%
  annotate_figure(
    left = text_grob("Number of days when the animal\nis detectable before being infectious",  rot = 90),
    bottom = text_grob("Sensitivity of diagnostic test for asymptomatics")
  )
ggsave(paste0("figures/final_param_set/fig3_p95.png"), bg="white")


a1 <- ggplot(results %>% filter(metric=="duration"), aes(x = se, y = days_detectable, fill = p95)) +
  geom_tile() +
  scale_fill_viridis_c() +
  # geom_text(aes(label = round(avg, 1)), size = 3) +
  facet_wrap(~effectiveness_insecticide) +
  theme_bw() +
  labs(
    x = "Sensitivity of diagnostic test for asymptomatics",
    y = "Number of days when the animal\nis detectable before being infectious",
    fill = "Value",
    title="Time to virus extinction (95th percentile)"
  )

a2 <- ggplot(results %>% filter(metric=="Ic"), aes(x = se, y = days_detectable, fill = p95)) +
  geom_tile() +
  scale_fill_viridis_c() +
  # geom_text(aes(label = round(avg, 1)), size = 3) +
  facet_wrap(~effectiveness_insecticide) +
  theme_bw() +
  labs(
    x = "Sensitivity of diagnostic test for asymptomatics",
    y = "Number of days when the animal\nis detectable before being infectious",
    fill = "Value",
    title="Total number of infected cows (95th percentile)"
  )

ggarrange(
  a1 + theme(axis.title= element_blank()),
  a2 + theme(axis.title = element_blank()),
  nrow = 2, ncol = 1,
  align = "hv", labels=c('A','B')) %>%
  annotate_figure(
    left = text_grob("Number of days when the animal\nis detectable before being infectious",  rot = 90),
    bottom = text_grob("Sensitivity of diagnostic test for asymptomatics")
  )
ggsave("figures/final_param_set/fig4_p95.png", bg="white")
