#-------------------------#
#--- analise MODELPERF ---#
#-------------------------#
#--- contact: clarabazzo@uni-bonn.de
library(stringr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(stringr)
library(multcompView)
library(AICcmodavg)


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

#--- get average across repeteitions
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
alldata_mean = cbind(alldata_mean, aggregate(reformulate(l_agg, 'rrmse'), alldata, mean)[,ncol(alldata_mean)])
colnames(alldata_mean)[length(colnames(alldata_mean))] = 'rrmse'

#--- here's the mean r2/rmse
alldata_mean


#-----------------------------#
#--- now read rank results ---#
#-----------------------------#

#--- remove alldata to relase memory
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


#---Select data for 'subset' = test, window & direction, 
#Para tarf_var = 'biomass' usar Window = '15x15' e direction '0'
#Para tarf_var = 'N_Species' usar Window = '7x7' e direction '135'
biomassdata <- alldata[alldata$targ_var == 'N_species'&
                         alldata$window_size == '7x7' &
                         alldata$direction == '135' &
                         alldata$subset == 'test'& 
                          alldata$data_filter == 'alldata'&
                         alldata$data_filter != 'alldata.csv',]

#--- make pretty labels
biomassdata$Classes = biomassdata$class_features
biomassdata$Model    = biomassdata$regression_model
biomassdata$Value    = biomassdata$r2
biomassdata$variable  = biomassdata$targ_var

#--- get unique lists
l_Treatment    = unique(biomassdata$Treatment)
l_Classes    = unique(biomassdata$Classes)
l_Model   = unique(biomassdata$Model)
l_Var     = unique(biomassdata$variable)

# Visualizando os valores únicos na coluna
valores_unicos <- unique(allrank$Var.Names)

# Imprimindo os valores únicos
print(valores_unicos)

#--- functions
ow_anova_tukey = function(df_test,
                          targ_var,
                          by_var){
  
  #--- perform an one-way ANOVA+Tukey Test
  #--- sources: 
  #---    - https://www.scribbr.com/statistics/anova-in-r/
  #---    - https://r-graph-gallery.com/84-tukey-test.html
  
  
  df_test$targ_var = df_test[,targ_var]
  df_test$by_var   = df_test[,by_var]
  
  #--- get ANOVA
  ow_anova = aov(targ_var ~ by_var, data = df_test)
  
  #--- get Tukey
  tuk = TukeyHSD(ow_anova)
  
  #--- get group letters
  if(length(tuk[['by_var']][,4]) >1){
    
    # Extract labels and factor levels from Tukey post-hoc 
    Tukey.levels <- tuk[['by_var']][,4]
    Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
    
    # I need to put the labels in the same order as in the boxplot :
    Tukey.labels$treatment=rownames(Tukey.labels)
    Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
    row.names(Tukey.labels) = NULL
    
  }else{
    
    #--- if only two treatments, either a or b
    
    Tukey.labels = 
      data.frame(Letters = 'a',
                 treatment = unique(df_test$by_var))
    
    #--- check for statistically difference
    if(tuk[['by_var']][,4] < 0.05){
      Tukey.labels$Letters[2]='b'
    }
    
  }
  
  return(list(tukey  = tuk,
              labels = Tukey.labels))
}


#-----------------------------#
#--- Run statistical tests ---#
#-----------------------------#

tuk_results = list()
for(v in l_Var){
  
  #--- filter data
  df_test = biomassdata[biomassdata$variable == v,]
  
  #--- target col
  targ_var = 'Value'
  
  #--- treatment
  by_var   = 'Model'
  
  tuk_labels = list()
  for(Treatment in l_Treatment){
    for(Classes in l_Classes){
      
      #--- filter except Model
      f = df_test$Treatment == Treatment & df_test$Classes == Classes
      
      #--- get tukey results
      tuk_res = ow_anova_tukey(df_test[f,], targ_var, by_var)
      
      #--- get labels and tag them
      tuk_lab = tuk_res$labels
      tuk_lab$Model = tuk_lab$treatment
      tuk_lab$Treatment  = Treatment
      tuk_lab$Classes  = Classes
      tuk_lab$treatment = by_var
      
      tuk_labels[[length(tuk_labels)+1]] = tuk_lab
      
    }
  }
  tuk_labels_Model = do.call(rbind,tuk_labels)
  
  #--- treatment
  by_var   = 'Classes'
  
  tuk_labels = list()
  for(Treatment in l_Treatment){
    for(Model in l_Model){
      
      #--- filter except Classes
      f = df_test$Treatment == Treatment & df_test$Model == Model
      
      #--- get tukey
      tuk_res = ow_anova_tukey(df_test[f,], targ_var, by_var)
      
      #--- get labels and tag them
      tuk_lab = tuk_res$labels
      tuk_lab$Classes = tuk_lab$treatment
      tuk_lab$Treatment  = Treatment
      tuk_lab$Model  = Model
      tuk_lab$treatment = by_var
      
      tuk_labels[[length(tuk_labels)+1]] = tuk_lab
      
    }
  }
  tuk_labels_Classes = do.call(rbind,tuk_labels)
  
  #--- merge both
  tuk_labels = 
    merge(tuk_labels_Classes,
          tuk_labels_Model,
          by = c('Treatment','Model','Classes'), 
          suffixes = c('_Classes','_Model'))
  
  #--- tag variable
  tuk_labels$variable = v
  
  #--- append
  tuk_results[[length(tuk_results)+1]] = tuk_labels
  
}
tuk_results = do.call(rbind, tuk_results)


#--- function to calculate the upper whisker
UW = function(x){
  #--- get upper whisker to place the letter
  Q1 = quantile(x, 0.25)
  Q3 = quantile(x, 0.75)
  
  upper_whisker = Q3 + (Q3 - Q1) * 1.5
  return(upper_whisker)
}

#--- calculate the upper whiskers for each combination of variables
biomassdata_UW = aggregate(Value ~ Treatment + Classes + Model + variable, biomassdata, UW)

#--- merge Tukey results with upper whisker data
tuk_letters = 
  merge(tuk_results,
        biomassdata_UW, 
        by = c('Treatment','Classes','Model','variable'))

#--- combine Tukey letters for each treatment
#--- assuming tuk_letters$Letters is already present in tuk_results
tuk_letters$Letters = paste0(toupper(tuk_letters$Letters_Classes),'',tuk_letters$Letters_Model)

#--------------------------------#
#--- Organize data and Plot ---#
#--------------------------------#

#--- font size
letter_size = 3

# Set dodge width according to the width of the boxplots
dodge_width <- 0.60


# Rename Treatment factor levels with new names
biomassdata$Treatment <- factor(biomassdata$Treatment,
                                levels = c("1", "2", "3", 'all'),
                                labels = c("Two Cut System", "Three Cut System", "Four Cut System", "Pooled Data"))

# Rename Classes factor levels with new names
biomassdata$Classes <- factor(biomassdata$Classes,
                              levels = c('dem', 'dem.glcm', 'dem.vi', 'dem.vi.glcm', 'glcm', 'vi', 'vi.glcm'),
                              labels = c('CH', 'CH.GLCM', 'CH.VI', 'CH.VI.GLCM', 'GLCM', 'VI', 'VI.GLCM'))

# Update tuk_letters with new names for Classes and Treatment
tuk_letters$Treatment <- factor(tuk_letters$Treatment,
                                levels = c("1", "2", "3", 'all'),
                                labels = c("Two Cut System", "Three Cut System", "Four Cut System", "Pooled Data"))

tuk_letters$Classes <- factor(tuk_letters$Classes,
                              levels = c('dem', 'dem.glcm', 'dem.vi', 'dem.vi.glcm', 'glcm', 'vi', 'vi.glcm'),
                              labels = c('CH', 'CH.GLCM', 'CH.VI', 'CH.VI.GLCM', 'GLCM', 'VI', 'VI.GLCM'))

# Define the desired order for Classes
levels_classes <- c('CH', 'VI', 'GLCM', 'CH.VI', 'CH.GLCM', 'VI.GLCM', 'CH.VI.GLCM')

# Update Classes variable in the dataframe to reflect the new order
biomassdata$Classes <- factor(biomassdata$Classes, levels = levels_classes)

# Ensure that the Classes variable in tuk_letters is also in the same order
tuk_letters$Classes <- factor(tuk_letters$Classes, levels = levels_classes)

# Define the desired order for Treatment
levels_treatment <- c("Two Cut System", "Three Cut System", "Four Cut System","Pooled Data")

# Update Treatment variable in the dataframe to reflect the new order
biomassdata$Treatment <- factor(biomassdata$Treatment, levels = levels_treatment)

# Ensure that the Treatment variable in tuk_letters is also in the same order
tuk_letters$Treatment <- factor(tuk_letters$Treatment, levels = levels_treatment)


#---- Create the ggplot

gg_boxplot <- ggplot(biomassdata, aes(y=Value*100, x=Model, group=interaction(Model,Classes), fill=Classes)) + 
  geom_boxplot(lwd=0.3, outlier.alpha = 1, outlier.color = NULL, outlier.size = 0.5, outlier.shape = 4) + 
  facet_grid(variable~Treatment, scales = 'free') + 
  theme_bw() + 
  theme(legend.position = 'bottom',
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text()) +
  geom_text(data=tuk_letters, aes(label=Letters, y=Value), 
            position=position_dodge(width=dodge_width), hjust=0.5, size=letter_size, vjust=-0.5) +
  scale_fill_brewer(palette = "Paired") + 
  labs(fill = "Classes",
       y = 'rRMSE(%)')

# Print the plot
print(gg_boxplot)


# Display the plot
print(gg_boxplot)

# Save the plot
ggsave(paste0(wd,'Results/boxplot.png'), gg_boxplot, dpi = 500, width = 12, height = 5)

#---- Create dotplot|boxplot

# Install and load the gghalves library
if (!require(gghalves)) install.packages("gghalves")
library(gghalves)

# Set the number of rows and columns for the facet grid
n_rows <- 2 # replace with the desired number of rows
n_cols <- 2 # replace with the desired number of columns

# Set a smaller dodge width for the points
dodge_width_points <- 0.85  # Adjust value as needed for points

# Set box width directly in geom_half_boxplot
box_width <- 0.75  # Adjust value as needed for boxes


# Create the plot

gg_half_box_dot <- ggplot(biomassdata, aes(x = Model, y = Value, fill = Classes)) + 
  facet_grid(facets = vars(variable, Treatment)) +
  geom_half_boxplot(aes(fill = Classes), position = position_dodge(width = dodge_width_points), outlier.shape = NA, width = box_width) + 
  geom_half_point(aes(colour = Classes), position = position_dodge(width = dodge_width_points)) + 
  theme_bw() +
  theme(legend.position = 'right',
        legend.text = element_text(size = 12),  # Increase font size in the legend
        legend.title = element_text(size = 14, face = "bold"),  # Increase size and make the legend title bold
        legend.key.size = unit(1.5, "lines"),  # Increase size of legend items
        axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y.left = element_text(size = 12, face = "bold"),  # Adjust left strip text size and font
        strip.text.y.right = element_text(size = 12, face = "bold"),  # Adjust right strip text size and font
        strip.placement = "outside") +
  geom_text(data = tuk_letters, aes(label = Letters, y = Value), 
            position = position_dodge(width = dodge_width), hjust = 0.5, size = letter_size, vjust = -0.5) +
  scale_fill_brewer(palette = "Paired") +  
  scale_colour_brewer(palette = "Paired") + 
  labs(fill = "Classes",
       colour = "Classes", 
       y = expression(R^2)) +
  ylim(0, 1) # Set y-axis limits

# View the plot
print(gg_half_box_dot)

# Save the plot
ggsave(filename = paste0(wd, 'Results/Nspecies_R2.png'), 
       plot = gg_half_box_dot, 
       dpi = 600, 
       width = 12, 
       height = 10)

