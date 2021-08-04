# Script for analysing b-catenin data obtained with R package cellPixels +++
# Author: Kai Budde
# Created: 2021/06/22
# Last changed: 2021/08/04

rm(list=ls())

# Input parameters #########################################################

data_directory <- "data"
plot_directory <- "plots"
output_directory <- "bcat"
file <- "bcat_intensity_summary.csv"
quality_result_to_be_excluded <- c("bad")

# General function parameters ##############################################
options(stringsAsFactors = FALSE, warn=-1)

library(dplyr)
library(ggplot2)

# Data input, cleaning, and additional columns #############################

input_file <- file.path(data_directory, file)
df_data <- read.csv(file = input_file, header = T)
output_dir <- paste(plot_directory, output_directory, sep="/")
dir.create(output_dir, showWarnings = FALSE)

# Load df_results
#load(paste(directory,"df_results.Rda", sep="/"))

# Delete all rows with specific description in quality check column
df_data <- df_data[(!df_data$manual_quality_check %in% quality_result_to_be_excluded),]

# Add column with information of magnification of the objective of the microscope
df_data$magnification <- NA
df_data$magnification <- gsub(pattern = ".+_([[:digit:]]+x)_.+",
                              replacement = "\\1",
                              x = df_data$fileName)

# Add column with information about which group an experiment belongs to
df_data$experiment_number <- NA
df_data$experiment_number <- gsub(pattern = "^V([[:digit:]]+)_.+",
                                  replacement = "\\1",
                                  x = df_data$fileName)

# Add column with information about which group an experiment belongs to
df_data$experiment_group <- NA
df_data$experiment_group[
  grepl(pattern = "stim", x = df_data$fileName, ignore.case = TRUE)] <- "Stimulation"
df_data$experiment_group[
  grepl(pattern = "contr", x = df_data$fileName, ignore.case = TRUE)] <- "Control"

# Add column with information about which position relative to the
# electrodes an image belongs to
df_data$image_position <- NA
df_data$image_position <- gsub(pattern = ".+-([MCP])[[:digit:]]+\\.czi",
                               replacement = "\\1",
                               x = df_data$fileName)

# Add column with information about which position relative to the
# electrodes an image belongs to
df_data$end_time <- NA
df_data$end_time <- gsub(pattern = ".+_([[:digit:]]+min).+",
                               replacement = "\\1",
                               x = df_data$fileName)
df_data$end_time <- factor(df_data$end_time, levels = c("30min", "60min", "120min"))


# Calculations of additional values to be added to the df ##################

# Beta-catenin: corrected total fluorescence inside the nucleus and inside
# the entire cell
df_data$bcat_corrected_total_fluorescence_nucleus <- NA
df_data$bcat_corrected_total_fluorescence_nucleus <-
  df_data$intensity_sum_red_nucleus_region -
  (df_data$intensity_mean_red_background*df_data$number_of_pixels_nucleus_region)

df_data$bcat_corrected_total_fluorescence_nucleus_per_no_of_nuclei <- NA
df_data$bcat_corrected_total_fluorescence_nucleus_per_no_of_nuclei <-
  df_data$bcat_corrected_total_fluorescence_nucleus /
  df_data$number_of_nuclei


df_data$bcat_corrected_total_fluorescence_cell <- NA
df_data$bcat_corrected_total_fluorescence_cell <-
  df_data$intensity_sum_red_foreground -
  (df_data$intensity_mean_red_background*df_data$number_of_pixels_foreground)

df_data$bcat_corrected_total_fluorescence_cell_per_no_of_nuclei <- NA
df_data$bcat_corrected_total_fluorescence_cell_per_no_of_nuclei <-
  df_data$bcat_corrected_total_fluorescence_cell /
  df_data$number_of_nuclei


# Beta-actin: corrected total fluorescence inside the entire cell
df_data$bactin_corrected_total_fluorescence_cell <- NA
df_data$bactin_corrected_total_fluorescence_cell <-
  df_data$intensity_sum_green_foreground -
  (df_data$intensity_mean_green_background*df_data$number_of_pixels_foreground)


df_data$bactin_corrected_total_fluorescence_cell_per_no_of_nuclei <- NA
df_data$bactin_corrected_total_fluorescence_cell_per_no_of_nuclei <-
  df_data$bactin_corrected_total_fluorescence_cell /
  df_data$number_of_nuclei


# df_data$intensity_sum_red_nucleus_region_per_no_of_nuclei <- NA
# df_data$intensity_sum_red_nucleus_region_per_no_of_nuclei <-
#   df_data$intensity_sum_red_nucleus_region / df_data$number_of_nuclei
#
# df_data$intensity_sum_red_outside_of_nucleus_region_per_no_of_nuclei <- NA
# df_data$intensity_sum_red_outside_of_nucleus_region_per_no_of_nuclei <-
#   df_data$intensity_sum_red_without_nucleus_region / df_data$number_of_nuclei

# Beta-actin
# df_data$intensity_sum_green_nucleus_region_per_no_of_nuclei <- NA
# df_data$intensity_sum_green_nucleus_region_per_no_of_nuclei <-
#   df_data$intensity_sum_green_nucleus_region / df_data$number_of_nuclei
#
# df_data$intensity_sum_green_outside_of_nucleus_region_per_no_of_nuclei <- NA
# df_data$intensity_sum_green_outside_of_nucleus_region_per_no_of_nuclei <-
#   df_data$intensity_sum_green_without_nucleus_region / df_data$number_of_nuclei
#
# df_data$intensity_sum_green_per_no_of_nuclei <- NA
# df_data$intensity_sum_green_per_no_of_nuclei <-
#   df_data$intensity_sum_green_full / df_data$number_of_nuclei


# Grouping results #########################################################

# Beta-Catenin inside nucleus
print("Beta-Catenin above/in nuclei per number of nuclei")
df_red_nuc <- group_by(df_data, magnification, end_time, experiment_number, experiment_group) %>%
  summarise(mean_intensity = mean(bcat_corrected_total_fluorescence_nucleus_per_no_of_nuclei))

df_red_nuc_2 <- df_red_nuc %>% tidyr::spread(key = experiment_group, value = mean_intensity)
df_red_nuc_2$stim_over_control <- df_red_nuc_2$Stimulation / df_red_nuc_2$Control

# Beta-Catenin inside entire cells
print("Beta-Catenin inside the entire cell per number of nuclei")
df_red_cell <- group_by(df_data, magnification, end_time, experiment_number, experiment_group) %>%
  summarise(mean_intensity = mean(bcat_corrected_total_fluorescence_cell_per_no_of_nuclei))

df_red_cell_2 <- df_red_cell %>% tidyr::spread(key = experiment_group, value = mean_intensity)
df_red_cell_2$stim_over_control <- df_red_cell_2$Stimulation / df_red_cell_2$Control


# Beta-Actin inside entire cells
print("Beta-Actin inside the entire cell per number of nuclei")
df_green_cell <- group_by(df_data, magnification, end_time, experiment_number, experiment_group) %>%
  summarise(mean_intensity = mean(bactin_corrected_total_fluorescence_cell_per_no_of_nuclei))

df_green_cell_2 <- df_green_cell %>% tidyr::spread(key = experiment_group, value = mean_intensity)
df_green_cell_2$stim_over_control <- df_green_cell_2$Stimulation / df_green_cell_2$Control


# Beta-Actin
# print("Beta-Actin inside/above nuclei per number of nuclei")
# df_green_nuc <- group_by(df_data, magnification, experiment_group) %>%
#   summarise(mean_intensity = mean(intensity_sum_green_nucleus_region_per_no_of_nuclei))
# print(df_green_nuc)
#
# print("Beta-Actin outside of nuclei per number of nuclei")
# df_green_cyt <- group_by(df_data, magnification, experiment_group) %>%
#   summarise(mean_intensity = mean(intensity_sum_green_outside_of_nucleus_region_per_no_of_nuclei))
# print(df_green_cyt)
#
# print("Beta-Actin in entire image per number of nuclei")
# df_green <- group_by(df_data, magnification, experiment_group) %>%
#   summarise(mean_intensity = mean(intensity_sum_green_per_no_of_nuclei))
# print(df_green)


# Plotting #################################################################

# Number (no) of nuclei per image (blue)
plot_no_nuclei <- ggplot(df_data, aes(x=end_time,
                                      y=number_of_nuclei,
                                      fill=experiment_group)) +
  geom_boxplot() +
  #facet_wrap(~end_time, scale="free") +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(size=18),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Number of nuclei per image") +
  xlab("") +
  scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
#print(plot_no_nuclei)

ggsave(filename = paste(output_dir, "number_of_nuclei.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir, "number_of_nuclei.png", sep="/"),
       width = 297, height = 210, units = "mm")



# Beta-Catenin (red)

# Beta-Catenin inside/above nuclei in entire image per number of nuclei
plot_red_nuc <- ggplot(df_red_nuc_2, aes(x=end_time,
                                      y=stim_over_control,
                                      color=experiment_number)) +
  geom_point(size = 6) +
  ylim(0.5,2) +
  theme_bw(base_size = 24) +
  theme(axis.title.y=element_text(size=18),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Ratio of corrected total fluorescence of Beta-Catenin\nabove nucleus per # of nuclei (Stim/Control)") +
  xlab("End time")

#print(plot_red_nuc)

ggsave(filename = paste(output_dir, "betacat_nuc.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir, "betacat_nuc.png", sep="/"),
       width = 297, height = 210, units = "mm")

# Beta-Catenin in cell in entire image per number of nuclei
plot_red_cell <- ggplot(df_red_cell_2, aes(x=end_time,
                                           y=stim_over_control,
                                           color=experiment_number)) +
  geom_point(size = 6) +
  ylim(0.5,2) +
  theme_bw(base_size = 24) +
  theme(axis.title.y=element_text(size=18),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Ratio of corrected total fluorescence of Beta-Catenin\nabove cell per # of nuclei (Stim/Control)") +
  xlab("End time")

#print(plot_red_cell)

ggsave(filename = paste(output_dir, "betacat_cell.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir, "betacat_cell.png", sep="/"),
       width = 297, height = 210, units = "mm")


# Beta-Actin in cell in entire image per number of nuclei
plot_green_cell <- ggplot(df_green_cell_2, aes(x=end_time,
                                           y=stim_over_control,
                                           color=experiment_number)) +
  geom_point(size = 6) +
  ylim(0.5,2) +
  theme_bw(base_size = 24) +
  theme(axis.title.y=element_text(size=18),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Ratio of corrected total fluorescence of Beta-Actin\nabove cell per # of nuclei (Stim/Control)") +
  xlab("End time")

#print(plot_green_cell)

ggsave(filename = paste(output_dir, "betaactin_cell.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir, "betaactin_cell.png", sep="/"),
       width = 297, height = 210, units = "mm")

# plot_red_nuc <- ggplot(df_data, aes(x=end_time,
#                                     y=intensity_sum_red_nucleus_region_per_no_of_nuclei,
#                                     fill=experiment_group)) +
#   geom_boxplot() +
#   #facet_wrap(~end_time, scale="free") +
#   theme_bw(base_size = 24) +
#   theme(axis.title.y=element_text(size=18),
#         #axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Sum of red fluorescence intensity above/\nin nucleus per number of nuclei per image") +
#   xlab("") +
#   scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
# print(plot_red_nuc)
#
# ggsave(filename = paste(output_dir, "betacat_nuc.pdf", sep="/"),
#        width = 297, height = 210, units = "mm")
# ggsave(filename = paste(output_dir, "betacat_nuc.png", sep="/"),
#        width = 297, height = 210, units = "mm")

#
# # Beta-Catenin inside/above nuclei in entire image per number of nuclei
# plot_red_cyt <- ggplot(df_data, aes(x=end_time,
#                                     y=intensity_sum_red_outside_of_nucleus_region_per_no_of_nuclei,
#                                     fill=experiment_group)) +
#   geom_boxplot() +
#   #facet_wrap(~end_time, scale="free") +
#   theme_bw(base_size = 24) +
#   theme(axis.title.y=element_text(size=18),
#         #axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Sum of red fluorescence intensity outside of\nnuclei regions per number of nuclei per image") +
#   xlab("") +
#   scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
# print(plot_red_cyt)
#
# ggsave(filename = paste(output_dir, "betacat_cyt.pdf", sep="/"),
#        width = 297, height = 210, units = "mm")
# ggsave(filename = paste(output_dir, "betacat_cyt.png", sep="/"),
#        width = 297, height = 210, units = "mm")
#
# # Beta-Actin (green)
#
# # Beta-Actin (green) outside of nuclei regions in entire image per number of nuclei
# plot_green_cyt <- ggplot(df_data, aes(x=end_time,
#                                       y=intensity_sum_green_outside_of_nucleus_region_per_no_of_nuclei ,
#                                       fill=experiment_group)) +
#   geom_boxplot() +
#   #facet_wrap(~magnification, scale="free") +
#   theme_bw(base_size = 24) +
#   theme(axis.title.y=element_text(size=18),
#         #axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Sum of green fluorescence intensity outside of\nnuclei regions per number of nuclei per image") +
#   xlab("") +
#   scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
# print(plot_green_cyt)
#
# ggsave(filename = paste(output_dir, "bactin_cyt.pdf", sep="/"),
#        width = 297, height = 210, units = "mm")
# ggsave(filename = paste(output_dir, "bactin_cyt.png", sep="/"),
#        width = 297, height = 210, units = "mm")
#
#
# # Beta-Actin (green) insinde of/above nuclei regions in entire image per number of nuclei
# plot_green_nuc <- ggplot(df_data, aes(x=end_time,
#                                       y=intensity_sum_green_nucleus_region_per_no_of_nuclei ,
#                                       fill=experiment_group)) +
#   geom_boxplot() +
#   #facet_wrap(~magnification, scale="free") +
#   theme_bw(base_size = 24) +
#   theme(axis.title.y=element_text(size=18),
#         #axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Sum of green fluorescence intensity above/\nin nucleus per number of nuclei per image") +
#   xlab("") +
#   scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
# print(plot_green_nuc)
#
# ggsave(filename = paste(output_dir, "bactin_nuc.pdf", sep="/"),
#        width = 297, height = 210, units = "mm")
# ggsave(filename = paste(output_dir, "bactin_nuc.png", sep="/"),
#        width = 297, height = 210, units = "mm")
#
#
# # Beta-Actin (green) in entire image per number of nuclei
# plot_green <- ggplot(df_data, aes(x=end_time,
#                                       y=intensity_sum_green_per_no_of_nuclei ,
#                                       fill=experiment_group)) +
#   geom_boxplot() +
#   #facet_wrap(~magnification, scale="free") +
#   theme_bw(base_size = 24) +
#   theme(axis.title.y=element_text(size=18),
#         #axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Sum of green fluorescence intensity\nin entire image per number of nuclei per image") +
#   xlab("") +
#   scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
# print(plot_green)
#
# ggsave(filename = paste(output_dir, "bactin_nuc.pdf", sep="/"),
#        width = 297, height = 210, units = "mm")
# ggsave(filename = paste(output_dir, "bactin_nuc.png", sep="/"),
#        width = 297, height = 210, units = "mm")

