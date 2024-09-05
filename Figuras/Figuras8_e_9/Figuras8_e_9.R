#################################################
###Author: Cássio Rampinelli & Saulo Aires
##################ANA-COMUC######################

#Clean R memory
rm(list=ls(all=TRUE))



############################################################################
###########################FIRST ATTEMPT TO GEV_MML#########################
############################################################################

############################################################################
####################################FIGURA 8###########################
############################################################################

#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras8_e_9")

folder_path="C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\IDF_Project\\Rio_Grande_do_Sul\\enviados_Cassio\\Figuras\\IDW_SA\\RMSE"

# List all .txt files in the folder
file_list <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)

# Initialize an empty list to store data frames
data_list <- list()

p<-c(0,0.15,0.25,0.3,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,9,10,15,20,30,40)
Tr_values<-c(100,2,25,5,500)
i=1
# Loop through each file and read the data
for (file in file_list) {
  # Read the data from the file
  data <- read.table(file, header = TRUE, sep = ",") # Adjust 'header' and 'sep' as needed
  
  Tr<-rep(Tr_values[i],length(p))
  data<-rbind(data,p)
  data<-rbind(data,Tr)
  
  rownames(data)[3]="p"
  rownames(data)[4]="Tr"
  
  # Store the data frame in the list
  data_list[[file]] <- data
  
  i=i+1
}

final.data.frame<-t(data_list[[1]])

#Combining all data frames
 for(i in 2:length(data_list)){
   
   final.data.frame<-rbind(final.data.frame,t(data_list[[i]]))
   
 }


final.data.frame=as.data.frame(final.data.frame)
final.data.frame$Tr<-as.factor(final.data.frame$Tr)
str(final.data.frame)

library(ggplot2)
library(dplyr)
library(tibble)

# Define custom line types and point shapes for each Tr value
line_types <- c("2" = "solid", "5" = "dashed", "25" = "dotted", "100" = "dotdash", "500" = "twodash")
point_shapes <- c("2" = 16, "5" = 17, "25" = 18, "100" = 19, "500" = 15)

# Define nudge_y values manually
nudge_values <- tibble::tribble(
  ~Tr, ~nudge_y,
  "2", -0.15,
  "5", -0.35,
  "25", -0.05,
  "100", -0.15,
  "500", -0.15
)

# Join nudge values to the data
final.data.frame <- final.data.frame %>%
  left_join(nudge_values, by = "Tr")

# Create the base plot
plot <- ggplot(final.data.frame, aes(x = p, y = RMSE, color = Tr, linetype = Tr, shape = Tr)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c("2" = "brown", "5" = "black", "25" = "darkgreen", "100" = "purple", "500" = "red")) +
  scale_linetype_manual(values = line_types) +
  scale_shape_manual(values = point_shapes) +
  labs(title = "p X RMSE",
       x = "p",
       y = "RMSE",
       color = "Tempo de Retorno",
       linetype = "Tempo de Retorno",
       shape = "Tempo de Retorno") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.width = unit(1.5, "cm"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10))

# Prepare data for labels with custom nudge_y values
label_data <- final.data.frame %>%
  group_by(Tr) %>%
  slice(which.min(RMSE)) %>%
  ungroup()

# Add labels for minimum values with customized nudge_y
plot <- plot + geom_text(data = label_data,
                         aes(label = paste("Min:", round(RMSE, 2))),
                         size = 3,
                         vjust = 0.5,
                         nudge_y = label_data %>%
                           left_join(nudge_values, by = "Tr") %>%
                           pull(nudge_y))

print(plot)



ggsave(filename = paste0("Figure8.png"), plot = plot, width = 10, height = 8, dpi = 500)




############################################################################
####################################FIGURA 9###########################
############################################################################



library(ggplot2)
library(dplyr)
library(tibble)

library(ggplot2)
library(dplyr)
library(tibble)


# Define custom point shapes for each Tr value
point_shapes <- c("2" = 16, "5" = 17, "25" = 18, "100" = 19, "500" = 15)

# Define nudge_x and nudge_y values manually
nudge_values <- tibble::tribble(
  ~Tr, ~nudge_x, ~nudge_y,
  "2", -0.15, -0.15,
  "5", -0.15, -0.35,
  "25", -0.15, -0.05,
  "100", -0.15, 0.15,
  "500", -0.15, -0.15
)

# Join nudge values to the data
final.data.frame <- final.data.frame %>%
  left_join(nudge_values, by = "Tr")

# Create the base plot
plot <- ggplot(final.data.frame, aes(x = RMSE, y = `%Dif.Slope`, color = Tr, shape = Tr)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("2" = "brown", "5" = "black", "25" = "darkgreen", "100" = "purple", "500" = "red")) +
  scale_shape_manual(values = point_shapes) +
  labs(title = " RMSE X %ΔS ",
       x = "RMSE",
       y = expression("%" * Delta * S),
       color = "Tempo de Retorno",
       shape = "Tempo de Retorno") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.width = unit(1.5, "cm"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10))

# Prepare data for labels with custom nudge_x and nudge_y values
label_data <- final.data.frame %>%
  group_by(Tr) %>%
  slice(which.min(RMSE)) %>%
  ungroup()

# Add labels for the p value associated with the minimum RMSE with customized nudge_x and nudge_y
plot <- plot + geom_text(data = label_data,
                         aes(label = paste("p:", round(p, 2))),
                         size = 4.5,
                         vjust = 0.5,
                         nudge_x = label_data %>%
                           left_join(nudge_values, by = "Tr") %>%
                           pull(nudge_x),
                         nudge_y = label_data %>%
                           left_join(nudge_values, by = "Tr") %>%
                           pull(nudge_y))

print(plot)


                         



ggsave(filename = paste0("Figure9.png"), plot = plot, width = 10, height = 8, dpi = 500)













