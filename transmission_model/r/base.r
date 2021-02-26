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
library(posterior)
library(here)

source(here("transmission_model/utils/geom-stepribbon.r"))
source(here("transmission_model/utils/gammaAlt.r"))

# check if sivep data exists if not download and read else just read
RELEASE_DATE_i = "08-02-2021"
SIVEP_ORIGINAL20 <- here(paste0("transmission_model/data/INFLUD-",RELEASE_DATE_i,".csv"))
SIVEP_ORIGINAL21 <- here(paste0("transmission_model/data/INFLUD21-",RELEASE_DATE_i,".csv"))

if(!file.exists(SIVEP_ORIGINAL20)){
  download.file(paste0("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2020/INFLUD-",RELEASE_DATE_i,".csv"),SIVEP_ORIGINAL20)
}
if(!file.exists(SIVEP_ORIGINAL21)){
  download.file(paste0("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2021/INFLUD21-",RELEASE_DATE_i,".csv"),SIVEP_ORIGINAL21)
}

df_SIVEP_ORIGINAL20 <- suppressMessages(read_csv2(SIVEP_ORIGINAL20))
df_SIVEP_ORIGINAL21 <- suppressMessages(read_csv2(SIVEP_ORIGINAL21))

df_SIVEP_original <- rbind(df_SIVEP_ORIGINAL20,df_SIVEP_ORIGINAL21)
deaths_nowcast_original = read.csv(here("transmission_model/data/nowcasted_daily_Manaus_class4and5_2021-02-11.csv"),
                                    stringsAsFactors = FALSE)
df_pop_original = read.csv(here("transmission_model/data/brazil_population.csv"),sep=";", stringsAsFactors = FALSE)
serial.interval = read.csv(here("transmission_model/data/serial_interval.csv"))
onset_paras_original = read.csv(here("transmission_model/data/statelevel_OnsetDeathParams.csv"))
pcr_and_sero <-read.csv(here("transmission_model/data/pcr_and_sero_pos.csv"))
pcr_genome_fraction <- read.csv(here("transmission_model/data/pcr_genome_fraction.csv"))

back_date = "2021-02-07" 
p1back_date = "2021-02-07" 
end_date = "2021-02-07" 
selected_state <- "MANAUS"
deaths_source <- "Nowcast4and5"
job = "base"
name="non_P1_P1_model"
rayleigh_par = 310
rel_IFR1 = 1/100
T2_date = ymd("2020-11-06")
ITER = 1000
WARM = 500
CORES = 3
THIN = 3
SEED = 4444
DELTA = 0.95
TREE = 12
JOBID=1

pcr_genome_fraction$date <- dmy(pcr_genome_fraction$date)
pcr_genome_fraction$negative <- pcr_genome_fraction$number - pcr_genome_fraction$positive
pcr_genome_fraction <- pcr_genome_fraction[pcr_genome_fraction$date > T2_date,]
pcr_genome_fraction <- pcr_genome_fraction[pcr_genome_fraction$date <= p1back_date,]

onset_paras = onset_paras_original[c("state","mean","cv")]  %>% drop_na()

deaths_nowcast <- deaths_nowcast_original
deaths_nowcast$"DateRep" = ymd(as.Date(deaths_nowcast$"Date"))
deaths_nowcast = deaths_nowcast[c("DateRep","Deaths","Var")]
deaths_nowcast$region = as.factor(selected_state)
deaths_nowcast_back = deaths_nowcast[deaths_nowcast$Date <= ymd(back_date),]
deaths_nowcast = deaths_nowcast[deaths_nowcast$Date <= ymd(end_date),]

if (deaths_source == "Nowcast4and5")
{
    df_SIVEP_original -> df_SIVEP
    df_SIVEP = df_SIVEP[df_SIVEP$ID_MUNICIP == selected_state,]
    df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5) || which(df_SIVEP$CLASSI_FIN==4),]
    df_SIVEP = df_SIVEP[which(df_SIVEP$EVOLUCAO==2),]
    df_SIVEP$DT_EVOLUCA = dmy(df_SIVEP$DT_EVOLUCA)
    filter_date = head(sort(df_SIVEP[df_SIVEP$CLASSI_FIN==5,]$DT_EVOLUCA),1)
    df_SIVEP = df_SIVEP %>% filter(DT_EVOLUCA >= ymd(filter_date))
    df_SIVEP = df_SIVEP[,c("DT_EVOLUCA", "ID_MUNICIP")]
    dim(df_SIVEP)
    df_SIVEP$Deaths = 1
    df_SIVEP = aggregate(. ~ID_MUNICIP+DT_EVOLUCA, data=df_SIVEP, sum, na.rm=TRUE)
    colnames(df_SIVEP) = c("region","DateRep","Deaths")
    df_SIVEP = rbind(df_SIVEP[df_SIVEP$DateRep < ymd(head(deaths_nowcast[c("DateRep","Deaths","region")],1)$DateRep),], deaths_nowcast[c("DateRep","Deaths","region")])
    regions <- unique(df_SIVEP$region)
    aux_df_zeros_brazil = NULL
    for( jj in 1:length(regions)){
      aux_sub = subset(df_SIVEP,region==regions[jj])
      aux_sub = aux_sub[order(aux_sub$DateRep),]
      last_date = aux_sub$DateRep[1]
      zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(last_date-1),"days"),
                         region=regions[jj],
                         Deaths=0)
      zeroes$DateRep <- ymd(zeroes$DateRep)
      aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
    }
    df_SIVEP_deaths=rbind.data.frame(aux_df_zeros_brazil,df_SIVEP)
} else
{
    df_SIVEP_original -> df_SIVEP
    df_SIVEP = df_SIVEP[df_SIVEP$ID_MUNICIP == selected_state,]
    df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5) || which(df_SIVEP$CLASSI_FIN==4),]
    df_SIVEP = df_SIVEP[which(df_SIVEP$EVOLUCAO==2),]
    df_SIVEP$DT_EVOLUCA = dmy(df_SIVEP$DT_EVOLUCA)
    df_SIVEP = df_SIVEP[,c("DT_EVOLUCA", "ID_MUNICIP")] #Date, and Federative unit patient residence.
    df_SIVEP$Deaths = 1
    df_SIVEP = aggregate(. ~ID_MUNICIP+DT_EVOLUCA, data=df_SIVEP, sum, na.rm=TRUE)
    colnames(df_SIVEP) = c("region","DateRep","Deaths")
    regions <- unique(df_SIVEP$region)
    aux_df_zeros_brazil = NULL
    for( jj in 1:length(regions)){
      aux_sub = subset(df_SIVEP,region==regions[jj])
      aux_sub = aux_sub[order(aux_sub$DateRep),]
      last_date = aux_sub$DateRep[1]
      zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(last_date-1),"days"),
                         region=regions[jj],
                         Deaths=0)
      zeroes$DateRep <- ymd(zeroes$DateRep)
      aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
    }
    df_SIVEP_deaths=rbind.data.frame(aux_df_zeros_brazil,df_SIVEP)
}
df_SIVEP_original -> df_SIVEP
df_SIVEP = df_SIVEP[df_SIVEP$ID_MUNICIP == selected_state,]
df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5),] 
df_SIVEP = df_SIVEP[,c("DT_SIN_PRI", "ID_MUNICIP")]
df_SIVEP$DT_SIN_PRI = dmy(df_SIVEP$DT_SIN_PRI)
df_SIVEP = df_SIVEP[order(df_SIVEP$DT_SIN_PRI),] %>% filter(DT_SIN_PRI >= ymd('2020-02-26')) #painel index case
df_SIVEP$Cases = 1
df_SIVEP = aggregate(. ~ID_MUNICIP+DT_SIN_PRI, data=df_SIVEP, sum, na.rm=TRUE)
colnames(df_SIVEP) = c("region","DateRep","Cases")
regions <- unique(df_SIVEP$region)
aux_df_zeros_brazil = NULL
for( jj in 1:length(regions)){
  aux_sub = subset(df_SIVEP,region==regions[jj])
  aux_sub = aux_sub[order(aux_sub$DateRep),]
  last_date = aux_sub$DateRep[1]
  zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(last_date-1),"days"),
                     region=regions[jj],
                     Cases=0)
  zeroes$DateRep <- ymd(zeroes$DateRep)
  aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
}
df_SIVEP_cases=rbind.data.frame(aux_df_zeros_brazil,df_SIVEP)
regions <- unique(df_SIVEP_deaths$region)
aux_df_zeros_brazil = NULL
for( jj in 1:length(regions)){
  zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(tail(sort(unique(df_SIVEP_deaths$DateRep)),1)),"days"),
                     region=regions[jj],
                     AUX=0)
  zeroes$DateRep <- ymd(zeroes$DateRep)
  aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
}
aux_df_zeros_brazil$key = paste0(aux_df_zeros_brazil$DateRep,aux_df_zeros_brazil$region)
df_SIVEP_deaths$key = paste0(df_SIVEP_deaths$DateRep,df_SIVEP_deaths$region)
df_SIVEP_cases$key = paste0(df_SIVEP_cases$DateRep,df_SIVEP_cases$region)
df_SIVEP_deaths = merge(aux_df_zeros_brazil, df_SIVEP_deaths, by = "key", all=TRUE)
df_SIVEP_cases = merge(aux_df_zeros_brazil, df_SIVEP_cases, by = "key", all=TRUE)
df_SIVEP = merge(df_SIVEP_deaths,df_SIVEP_cases, by = "key")
df_SIVEP = df_SIVEP[,c("DateRep.y.x","region.y.x","Deaths","Cases")]
colnames(df_SIVEP) = c("DateRep","region","Deaths","Cases")
df_SIVEP[c("Cases")] = df_SIVEP[c("Cases")] %>% replace(is.na(.), 0)
df_SIVEP = df_SIVEP %>% group_by(region) %>% mutate(cumulativeCases = cumsum(Cases))
df_SIVEP = df_SIVEP %>% group_by(region) %>% mutate(cumulativeDeaths = cumsum(Deaths))

df_pop <- df_pop_original[c("region","population")]
df_pop <- df_pop[1:(nrow(df_pop)-1),]
df_SIVEP <- merge(x = df_SIVEP, y = df_pop, by = "region", all = TRUE)
df_SIVEP <- df_SIVEP[order(df_SIVEP$DateRep),]

df=df_SIVEP %>% filter(DateRep <= ymd(back_date))
START_TIME=ymd(as.Date("2020-03-03"))
END_TIME=ymd(as.Date(end_date))
RANGE_TIME=seq(START_TIME,END_TIME,by = '1 day')
StanModel = as.character(job)
models = job
writeDir = job
d<-df
countries = selected_state
N2 = length(RANGE_TIME)
dates = list()
reported_cases = list()
stan_data = list(M=length(countries),
		 N=NULL,
                 deaths=NULL,
		 deaths_combined=NULL,
		 f=NULL,
		 N0=6,
		 cases=NULL,
		 SI=NULL,
		 EpidemicStart = NULL,
                 pop = NULL,
                 par = NULL,
                 T2 = NULL, 
                 phylo_N = NULL, 
		 phylo_PSamples = NULL, 
		 phylo_NSamples = NULL,
                 sero_N = NULL, 
		 sero_sigma = NULL, 
		 sero_prev = NULL,
                 len_excess = NULL, 
		 excess_N = NULL, 
		 excess = NULL,
                 nowcast_sd_N = NULL, 
		 nowcast_sd = NULL, 
		 nowcast_sd_len = NULL
		 ) 
deaths_by_country = list()
deaths_by_country_combined = list()
mean1 = 5.1; cv1 = 0.86;
x1 = rgammaAlt(1e6,mean1,cv1)
aux.epidemicStart = NULL
Country = selected_state
mean2 = onset_paras[onset_paras$state == Country,]$mean
cv2 = onset_paras[onset_paras$state == Country,]$cv
x2 = rgammaAlt(1e6,mean2,cv2)
ecdf.saved = ecdf(x1+x2)
d1=d[d$region==Country,]
d1$DateRep = as.Date(d1$DateRep)
d1_pop = df_pop[df_pop$region==Country,]
d1 = d1[order(d1$DateRep),]
index = which(d1$Cases>0)[1]
index1 = which(cumsum(d1$Deaths)>=5)[1] # 5, 10
index2 = index1-26 #26
d1=d1[index2:nrow(d1),]

aux.epidemicStart = c(aux.epidemicStart,d1$DateRep[index1+1-index2])
stan_data$EpidemicStart = c(stan_data$EpidemicStart,index1+1-index2)
stan_data$EpidemicStart = as.array(stan_data$EpidemicStart)
stan_data$pop = c(stan_data$pop, d1_pop$population)
stan_data$pop = as.array(stan_data$pop)
stan_data$par = rayleigh_par
dates[[Country]] = d1$DateRep
N = length(d1$Cases)
forecast = N2 - N
if(forecast < 0) {
  N2 = N
  forecast = N2 - N
}
convolution1 = function(u) (rel_IFR1 * ecdf.saved(u))
f = rep(0,N2) 
f[1] = (convolution1(1.5) - convolution1(0))
for(i in 2:N2)
{
      f[i] = (convolution1(i+.5) - convolution1(i-.5))
}
cases = as.vector(as.numeric(d1$Cases))
deaths = as.vector(as.numeric(d1$Deaths))
stan_data$N = c(stan_data$N,N)
stan_data$f = cbind(stan_data$f,f)
stan_data$deaths = cbind(stan_data$deaths,deaths)
stan_data$cases = cbind(stan_data$cases,cases)
stan_data$N2 = N2
stan_data$x=1:N2
if(length(stan_data$N) == 1) {
  stan_data$N = as.array(stan_data$N)
}
stan_data$W <- ceiling(stan_data$N2/7)
stan_data$week_index <- matrix(1,stan_data$M,stan_data$N2)
    for(state.i in 1:stan_data$M) {
    stan_data$week_index[state.i,] <- rep(2:(stan_data$W+1),each=7)[1:stan_data$N2]
    last_ar_week = which(dates[[state.i]]==max(df$Date) -  7)
    stan_data$week_index[state.i,last_ar_week:ncol(stan_data$week_index)] <-
    stan_data$week_index[state.i,last_ar_week]
}
stan_data$AR_SD_MEAN = 0.2
stan_data$P <- 0
stan_data$SI = serial.interval$fit[1:N2]
stan_data$T2 = seq(1, length(dates[[Country]]))[dates[[Country]] %in%  T2_date]

if(length(stan_data$deaths < N2))
{
    stan_data$deaths = rbind(stan_data$deaths,as.matrix(rep(-1,N2-N)))
}

stan_data$CasesStart <- 20 #40

length_pcr_sero_data <- length(pcr_and_sero$p_PCR_positive)
padding <- stan_data$N2 - length_pcr_sero_data
stan_data$PCR_pos_prob <- as.matrix(c(pcr_and_sero$p_PCR_positive, rep(0, padding)))
stan_data$seroconv_cdf <- as.matrix(c(pcr_and_sero$cum_seropositive, rep(1, padding)))
stan_data$serorev_surv <- as.matrix(1 - pweibull(seq(1, N2), shape = 2.933, scale = 208))

stan_data$phylo_N = seq(1, length(dates[[Country]]))[dates[[Country]] %in%  pcr_genome_fraction$date]
stan_data$phylo_N_len = length(stan_data$phylo_N)
stan_data$phylo_PSamples = pcr_genome_fraction$positive
stan_data$phylo_NSamples = pcr_genome_fraction$negative

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m = stan_model(here(paste0('transmission_model/stan_models/',StanModel,'.stan')))

fit = sampling(m,
               data=stan_data,
               iter=ITER,
               warmup=WARM,
               chains=CORES,
               cores=CORES,
	       thin=THIN,
	       seed=SEED,
               control = list(adapt_delta = DELTA, max_treedepth = TREE))
out = rstan::extract(fit)
filename <- paste0(StanModel,'-',JOBID)

if(!dir.exists(here(paste0("transmission_model/figures/",job)))){
  dir.create(here(paste0("transmission_model/figures/",job)), recursive = TRUE)
}
if(!dir.exists(here(paste0("transmission_model/results/",job)))){
  dir.create(here(paste0("transmission_model/results/",job)), recursive = TRUE)
}
save.image(here(paste0("transmission_model/results/",writeDir,"/",name,".Rdata")))
