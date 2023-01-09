library(readxl)
map_names <- excel_sheets("data/HeatMap.xlsx")     

map_list <- lapply(map_names, function(x) {          # Read all sheets to list
  
  as.data.frame(read_excel("data/HeatMap.xlsx",  sheet = x)) 
  
} )


names(map_list) <- map_names 

source("src/packages.R")
library(scales)

#Mauve
#D3B1C2

#Lavender
#C197D2

#Black
#211522

#Orchid
#613659
#pal <- colorRampPalette(c("orange1","darkorange3","coral3", "firebrick1")) 
pal <- colorRampPalette(c("lightblue1","mediumblue","orangered2"))
#pal <- colorRampPalette(c("lightblue","steelblue2","blue3"))


p1 <- ggplot(map_list$NAFE %>% filter(Stage == 'DEVELOPMENTAL'), aes(casefold(structure,upper = TRUE), Gene)) +
  geom_tile(aes(fill = Expression), colour = "white") +
  ggtitle("NAFE") +
  scale_fill_gradientn(colours = pal(50)) +
  labs(fill='Dev Score') +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2,face = "bold"),
        axis.text.y = element_text(size = rel(0.8),face = "italic")) + 
  scale_x_discrete(labels = wrap_format(20))

p2 <- ggplot(map_list$DEE %>% filter(Stage == 'DEVELOPMENTAL'), aes(casefold(structure, upper = TRUE), Gene)) +
  geom_tile(aes(fill = Expression), colour = "white") +
  ggtitle("DEE") +
  scale_fill_gradientn(colours = pal(50)) +
  labs(fill='Dev Score') +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2, size = 7,face = "bold"),
        axis.text.y = element_text(size = rel(0.8),face = "italic")) + 
  scale_x_discrete(labels = wrap_format(20))

p3 <- ggplot(map_list$GGE %>% filter(Stage == 'DEVELOPMENTAL'), aes(casefold(structure,upper = TRUE), Gene)) +
  geom_tile(aes(fill = Expression), colour = "white") +
  ggtitle("GGE") +
  scale_fill_gradientn(colours = pal(50)) +
  labs(fill='Dev Score') +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2,face = "bold"),
        axis.text.y = element_text(size = rel(0.8),face = "italic")) + 
  scale_x_discrete(labels = wrap_format(20))



patchwork = p1 + p2 + p3

patchwork[[1]] = patchwork[[1]] + theme(legend.position = "none") + 
  xlab('') 

patchwork[[2]] = patchwork[[2]] + theme(legend.position = "none",
) + ylab('') + xlab('')

patchwork[[3]] = patchwork[[3]]  + ylab('') +
  xlab('')

patchwork 
png("graph/Fig2B_ABA_Heatmap.png", width = 42, height = 21, units = 'cm', res = 300)
patchwork
dev.off()


png("graph/Fig2A_DEE_Heatmap.png", width = 42, height = 21, units = 'cm', res = 300)
p2
dev.off()
