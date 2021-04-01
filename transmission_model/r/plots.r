library(rstan)
library(matrixStats)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(scales)
library(tidyverse)
library(dplyr)
library(abind)
library(xtable)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(cowplot)
library(openxlsx)
library(loo)
library(cmdstanr)
library(ks)
library(ggpubr)
library(here)
library(optparse)
library(ggExtra)
option_list <- list(
    make_option(c("--fileName"),action="store", default=here("transmission_model/results/base/non_P1_P1_model.Rdata"),help="File to be loaded [default \"%default\"]"),
    make_option(c("--beta"),action="store", default="0.2",help="Cross immunity factor [default \"%default\"]")
)

opt <- parse_args(OptionParser(option_list=option_list))

load(opt$fileName)
out <- rstan::extract(fit)
posterior <- as.matrix(fit)

# posteriors of interest
######################################################### 
posteriors_of_interest <- 
    bind_rows(
        tibble(parameter = "Transmissibility Increase(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$R_difference - 1,probs = c(0.25)),0), round(100*quantile(out$R_difference - 1,probs = c(0.75)),0))
        ),
        tibble(parameter = "Relative Risk Increase(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$RR - 1,probs = c(0.25)),0), round(100*quantile(out$RR - 1,probs = c(0.75)),0))
        ),
        tibble(parameter = "Cross Immunity(%)",
               `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(1 - out$cross,probs = c(0.25)),0), round(100*quantile( 1 - out$cross,probs = c(0.75)),0))
        )
    )
print(posteriors_of_interest)
grid.table(posteriors_of_interest)
############################################################
# cross immunity vs transmissibility
############################################################
x <- data.frame(R_difference = out$R_difference, cross =out$cross)
kd <- ks::kde(x, compute.cont=TRUE,H=matrix(c(0.01,0,0,0.0033),2))
contour_95 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["5%"])[[1]]))
contour_75 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["25%"])[[1]]))
contour_50 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["50%"])[[1]]))
contour_25 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["75%"])[[1]]))
contour_5 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                              z=estimate, levels=cont["95%"])[[1]]))
pA <- mcmc_scatter(posterior, pars = c("R_difference","cross"), size = 3.5, alpha = 0.1) +
    ylab("Cross-immunity") +
    xlab("Transmissibility increase") +
    scale_y_continuous(expand=c(0,0)) +
    theme_pubr(base_size = 26) +
    ggplot2::geom_path(aes(x, y), data=contour_95) +
    ggplot2::geom_path(aes(x, y), data=contour_75) +
    ggplot2::geom_path(aes(x, y), data=contour_50) +
    ggplot2::geom_path(aes(x, y), data=contour_25) +
    ggplot2::geom_path(aes(x, y), data=contour_5) 
pA <- ggMarginal(pA, R_difference, cross, type = c("density"),
                 margins = 'both',
                 size = 4,
                 colour = "#589e73",
                 fill = '#589e73',
                 alpha = 0.25,
                 xparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)),
                 yparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)))
pA
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PA_contour.pdf")),plot=pA)
############################################################
# cross immunity vs relative risk increase
############################################################
x <- data.frame(R_difference = out$RR, cross = out$cross)
kd <- ks::kde(x, compute.cont=TRUE,H=matrix(c(0.015,0,0,0.0075),2))
contour_95 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["5%"])[[1]]))
contour_75 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["25%"])[[1]]))
contour_50 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["50%"])[[1]]))
contour_25 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                               z=estimate, levels=cont["75%"])[[1]]))
contour_5 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                              z=estimate, levels=cont["95%"])[[1]]))

pB <-  mcmc_scatter(posterior, pars = c("RR[1]","cross"), size = 3.5, alpha = 0.1) +
    ylab("Cross-immunity") +
    xlab("Relative risk of mortality") +
    scale_y_continuous(expand=c(0,0)) +
    theme_pubr(base_size = 26) +
    xlim(c(0,4)) +
    ggplot2::geom_path(aes(x, y), data=contour_95) +
    ggplot2::geom_path(aes(x, y), data=contour_75) +
    ggplot2::geom_path(aes(x, y), data=contour_50) +
    ggplot2::geom_path(aes(x, y), data=contour_25) +
    ggplot2::geom_path(aes(x, y), data=contour_5) 
pB <- ggMarginal(pB, R_difference, cross, type = c("density"),
                 margins = 'both',
                 size = 4,
                 colour = "#589e73",
                 fill = '#589e73',
                 alpha = 0.25,
                 xparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)),
                 yparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)))

pB
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PB_contour.pdf")),plot=pB)
############################################################
# Deaths Curve
############################################################
plt_deaths <- function(outcheck)
{
    excess= read.csv(here("transmission_model/data/excess_burials.csv"))
    colnames(excess) = c("Date","Deaths_raw","Mean7_Deaths","X")
    excess = excess[c("Date","Deaths_raw")]
    excess$Date = dmy(excess$Date)
    excess$Deaths_raw = as.numeric(excess$Deaths_raw)
    df_SIVEP_original -> df_SIVEP
    df_SIVEP = df_SIVEP[df_SIVEP$ID_MUNICIP == selected_state,]
    df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5) || which(df_SIVEP$CLASSI_FIN==4),]
    df_SIVEP = df_SIVEP[which(df_SIVEP$EVOLUCAO==2),]
    df_SIVEP$DT_EVOLUCA = dmy(df_SIVEP$DT_EVOLUCA)
    filter_date = head(sort(df_SIVEP[df_SIVEP$CLASSI_FIN==5,]$DT_EVOLUCA),1)#only after first confirmed
    df_SIVEP = df_SIVEP %>% filter(DT_EVOLUCA >= ymd(filter_date))
    df_SIVEP = df_SIVEP[,c("DT_EVOLUCA", "ID_MUNICIP")]
    dim(df_SIVEP)
    df_SIVEP$Deaths = 1
    df_SIVEP = aggregate(. ~ID_MUNICIP+DT_EVOLUCA, data=df_SIVEP, sum, na.rm=TRUE)
    sum(df_SIVEP$Deaths)
    colnames(df_SIVEP) = c("region","DateRep","Deaths")
    
    deaths_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.975)
    )
    
    deaths_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.025),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.975)
    )
    
    deaths_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v1[,,], probs=.75)
    )
    
    deaths_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "deaths" = colMeans(outcheck$E_deaths),
        "deaths" = colMeans(outcheck$E_deaths_v1),
        "deaths_l" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.25),
        "deaths_u" = colQuantiles(outcheck$E_deaths_v2[,,], probs=.75)
    )
    
    deaths_df = rbind(deaths_df95v1,deaths_df95v2,deaths_df50v1,deaths_df50v2)
    
    
    ggplot() + 
        
        geom_ribbon(data = deaths_df, 
                    aes(x=time, ymax=deaths_l, ymin=deaths_u, fill = key)) +
        scale_fill_manual(name = "", labels = c("50%","50%","95%","95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.45),
                                     alpha("tan2",0.45))) +
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 1)) +
        xlab("") +
        ylab("Daily deaths")+
        ylim(-5,160)+
        geom_point(data = excess, aes(x=Date, y = Deaths_raw, alpha=0.25))+
        geom_point(data = df_SIVEP, aes(x=DateRep, y=Deaths),col="firebrick2",alpha=0.5) +
        scale_x_date(date_breaks = "1 month", labels = date_format("%b")) +
        theme(axis.text.x = element_text(angle = 50, hjust = 1))  +
        theme(legend.position="none")
}

p_C <- plt_deaths(out)
p_C
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PC_contour.pdf")),plot=p_C)
############################################################
# P1 fraction curve
############################################################
plt_P1fraction <- function(outcheck)
{
    datapoints = readRDS(here("transmission_model/data/datapoints.rds"))
    colnames(datapoints)[1] <- "date"
    datapoints = datapoints[datapoints$date > ymd("2020-11-01"),]
    #     dates[[Country]] = c(dates[[Country]],seq(tail(dates[[Country]],1)+1,tail(dates[[Country]],1)+342-311,by = '1 day'))
    E_fraction_df_95 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("nintyfive", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.975),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.025)
    )
    E_fraction_df_50 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50", length(dates[[Country]])),
        "E_fraction" = colMeans(outcheck$E_fraction),
        "E_fraction_ui" = colQuantiles(outcheck$E_fraction[,,], probs=.75),
        "E_fraction_li" = colQuantiles(outcheck$E_fraction[,,], probs=.25)
    )
    E_fraction_df = rbind(E_fraction_df_95,E_fraction_df_50)
    late_E_fraction_df = E_fraction_df[E_fraction_df$time >= ymd("2020-11-06"),]
    
    ggplot() + 
        geom_point(data = datapoints, aes(x = date, y = prop)) +
        geom_errorbar(data = datapoints, aes(x = date, ymin=lower_CI, ymax=upper_CI,width=2))+
        geom_ribbon(data = late_E_fraction_df, aes(x=time, ymax=E_fraction_ui, ymin=E_fraction_li, fill = key))+
        
        scale_fill_manual(name = "", labels = c("50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("Deepskyblue4",0.55))) +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        xlab("") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab("P1 fraction") +
        
        
        scale_x_date(date_breaks = "month", labels = date_format("%e %b")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_D <- plt_P1fraction(out)
p_D
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PD_contour.pdf")),plot=p_D)
############################################################
# Seropervalence Curve
############################################################
plt_sero_conv <- function(outcheck)
{
    manausPopulation = 2219580
    seroconv_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("95v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.025)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.975)/stan_data$pop[[1]]
    )
    seroconv_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v1", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v1)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v1[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("50v2", length(dates[[Country]])),
        "seroconv" = colMeans(outcheck$seroconv_v2)/stan_data$pop[[1]],
        "seroconv_l" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.25)/stan_data$pop[[1]],
        "seroconv_u" = colQuantiles(outcheck$seroconv_v2[,,1], probs=.75)/stan_data$pop[[1]]
    )
    seroconv_df = rbind(seroconv_df95v1,seroconv_df95v2,seroconv_df50v1,seroconv_df50v2)
    
    
    manaus_seroprev = read.csv(here("transmission_model/data/manaus_seroprev.csv"))
    manaus_seroprev$date = ymd(manaus_seroprev$date)
    manaus_seroprev = manaus_seroprev[manaus_seroprev$upper < 150 & manaus_seroprev$upper > 0 ,] 
    manaus_seroprev$sero_sigma = (manaus_seroprev$prevalence - manaus_seroprev$lower)/2 
    tail(manaus_seroprev,1)
    
    p1=seroconv_df %>% ggplot() +
        geom_ribbon(aes(x = time, ymin = seroconv_l, ymax = seroconv_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.85),
                                     alpha("tan2",0.85),
                                     alpha("Deepskyblue4",0.55),
                                     alpha("tan2",0.55))) +
        geom_point(data=manaus_seroprev,aes(x=date,y=prevalence/100)) +
        geom_errorbar(data=manaus_seroprev,aes(x=date,ymin=lower/100, ymax=upper/100,width=10)) +
        
        xlab("") +
        ylab("Cumulative incidence per capita\n") +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.2, 0.9)) +
        scale_x_date(date_breaks = "1 month", labels = date_format("%b")) +
        theme(axis.text.x = element_text(angle = 50, hjust = 1))
    
    p1
}
p_E <-  plt_sero_conv(out)
p_E
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PE_contour.pdf")),plot=p_E)
############################################################
# Rt Curve
############################################################
plt_rt <- function(outcheck)
{
    rt_df95v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_95)", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.975)
    )
    rt_df95v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_95", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.025),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.975)
    )
    rt_df50v1 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v1_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v1),
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v1[,,], probs=.75)
    )
    rt_df50v2 = data.frame(
        "time" = dates[[Country]],
        "key" = rep("rt_v2_50", length(dates[[Country]])),
        "rt" = colMeans(outcheck$Rt_adj_immune_v2)/stan_data$pop[[1]],
        "rt_l" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.25),
        "rt_u" = colQuantiles(outcheck$Rt_adj_immune_v2[,,], probs=.75)
    )
    
    rt_df = rbind(rt_df50v1,rt_df50v2,rt_df95v1,rt_df95v2)
    tail(rt_df,1)
    rt_df %>% ggplot() +
        geom_ribbon(aes(x = time, ymin = rt_l, ymax = rt_u, fill =key)) +
        geom_ribbon(aes(x = time, ymin = rt_l, ymax = rt_u, fill =key)) +
        scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                          values = c(alpha("Deepskyblue4",0.55),
                                     alpha("Deepskyblue4",0.35),
                                     alpha("tan2",0.55),
                                     alpha("tan2",0.35))) +
        
        theme_pubr(base_size = 26) +
        theme(legend.position = c(0.25, 0.9)) +
        xlab("") +
        scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
        ylab(expression(R[t])) +
        
        scale_x_date(date_breaks = "1 month", labels = date_format("%b")) +
        theme(axis.text.x = element_text(angle = 50, hjust = 1)) 
}


p_F <- plt_rt(out)
p_F
ggsave(here(paste0("transmission_model/figures/",job,"/fig_PF_contour.pdf")),plot=p_F)
############################################################
library(cowplot)
plot_grid(plot_grid(pA, p_C, p_E, labels = c("A", "C", "E"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(1.05,1,1)), 
          plot_grid(pB, p_D, p_F, labels = c("B", "D", "F"),nrow = 1, align = "h",axis = "b",
                    rel_widths = c(1.05,1.1, 1)),
          nrow = 2 )
ggsave(here(paste0("transmission_model/figures/",job,"/fig_4.pdf")))
