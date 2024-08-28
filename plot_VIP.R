#---------------------------------#
#--- Plot Variable Importance ---#
#--------------------------------#
#--- contact: clarabazzo@uni-bonn.de
library(stringr)
library(ggplot2)
library(forcats)
library(dplyr)

wd = paste0(getwd(),'/')
MODELPERFzip = 'Results/ALLDATA/MODELPERF/MODELPERF_REPS_TREAT.zip'

#--- list files in zip
list_fn_inzip = unzip(paste0(wd, MODELPERFzip), list = T)
list_fn_inzip = list_fn_inzip[grep('MODELPERF_', list_fn_inzip$Name),]

#--- read all results inside .zip file
alldata = list()
for(fn in list_fn_inzip$Name){
  f = read.csv(unzip(paste0(wd, MODELPERFzip), fn), as.is = T)
  f$data_filter = str_split(fn, '_', simplify = T)[,4]
  alldata[[length(alldata)+1]] = f
}
alldata = do.call(rbind, alldata)

alldata$targ_var = NA
alldata$targ_var[grepl('N_species',alldata$model)] = 'N_species'
alldata$targ_var[grepl('biomass',alldata$model)] = 'biomass'

#--- get average across repetitions
l_agg = c('model',
          'n',
          'subset',
          'total_features',
          'regression_model',
          'class_features',
          'window_size',
          'direction',
          'data_filter',
          'Treatment',
          'targ_var')

alldata_mean = aggregate(reformulate(l_agg, 'r2'), alldata, mean)
alldata_mean = cbind(alldata_mean, aggregate(reformulate(l_agg, 'rmse'), alldata, mean)[,ncol(alldata_mean)])
colnames(alldata_mean)[length(colnames(alldata_mean))] = 'rmse'


#-----------------------------#
#------- Read rank results ---#
#-----------------------------#

#--- remove alldata to release memory
rm(alldata)
gc()

list_fn_inzip = unzip(paste0(wd, MODELPERFzip), list = T)
list_fn_inzip = list_fn_inzip[grep('MODELRANK_', list_fn_inzip$Name),]

#--- read all results inside .zip file
allrank = list()
for(fn in list_fn_inzip$Name){
  f = read.csv(unzip(paste0(wd, MODELPERFzip), fn), as.is = T)
  f$data_filter = str_split(fn, '_', simplify = T)[,4]
  f$targ_var = str_split(fn, '_', simplify = T)[,5]
  allrank[[length(allrank)+1]] = f
}
allrank = do.call(rbind, allrank)
gc()

allrank$targ_var[allrank$targ_var == 'N'] = 'N_species'


#--- Select data for window & direction, and targ_var
# For targ_var = 'biomass' use Window = '15x15' and direction '135'
# For targ_var = 'N_species' use Window = '7x7' and direction '135'
nspecies <- allrank[allrank$targ_var == 'N_species'&
                      allrank$window_size == '15x15' &
                      allrank$direction == '135' & 
                      allrank$data_filter == 'alldata'&
                      allrank$regression_model == 'RF'& #change the regression model do 'PLS' for PLS VIP
                      allrank$data_filter != 'alldata.csv',]

#--- Make pretty labels
nspecies$Var.Names <- sub("^dem", "CH", nspecies$Var.Names)
nspecies$Var.Names <- sub("^vi", "VI", nspecies$Var.Names)
nspecies$Var.Names <- sub("^glcm", "GLCM", nspecies$Var.Names)

#--- Prepare the data
nspecies_summary <- nspecies %>%
  group_by(Treatment, Var.Names) %>%
  summarise(Average_Relative_Importance = mean(relative_importance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Treatment = factor(Treatment, levels = c('1', '2', '3', 'all')))

#--- Filter the top 10 for each treatment and reorder variables according to relative importance
top10_by_treatment <- nspecies_summary %>%
  arrange(Treatment, desc(Average_Relative_Importance)) %>%
  group_by(Treatment) %>%
  slice_max(order_by = Average_Relative_Importance, n = 10) %>%
  ungroup() %>%
  mutate(Var.Names = fct_reorder(Var.Names, Average_Relative_Importance))

#--- Define better titles for treatments
treatment_titles <- c('1' = 'Two-Cut System', '2' = 'Three-Cut System', '3' = 'Four-Cut System', 'all' = 'Pooled Data')

#--- Create the plot
p <- ggplot(top10_by_treatment, aes(x = Var.Names, y = Average_Relative_Importance*10, fill = Treatment)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Treatment, scales = 'free_y', nrow = 1, labeller = as_labeller(treatment_titles)) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = 'Average Relative Importance', title = 'Top 10 Variables by Relative Importance for Number of Species across Treatments for RF') +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, angle = 0, hjust = 1, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, face = "bold"),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), # Adding bold to the plot title
        legend.text = element_text(face = "bold")) + 
  scale_fill_brewer(palette = "Paired") +  
  scale_colour_brewer(palette = "Paired")+
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(0, 100))

#--- Display the plot
print(p)

#--- Save the plot
ggsave("VIP_RF.png", path = paste0(wd,"Results/"), width = 16, height = 6, dpi = 300)


