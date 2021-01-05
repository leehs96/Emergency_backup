

path <- "/users/hslee/조선욱_교수님/201222/score/out/새 폴더/"
data2 <-list.files(path = path, pattern = '*.csv', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("_heat.csv","_scatter",b))
  a <- (read.csv(data2[i], header = T , sep = ',', row.names = 1))
  #a <- dplyr::select(a, grep('.01$', names(a)))
  a <- as.data.frame(t(a))
  assign(name, a)
}

setwd('~/TCGA_RNA_data')
normal <- t(read.csv('GTEX_TCGA_gene.csv', sep = ',', header = T, row.names = 1))
normal_ID <- read.csv('GTEx_Analysis_v8_Annotations_SampleAttributesDS.csv', sep = ',', header = T, row.names = 1 )

normal <- as.data.frame(normal)

normal_ID_CD14 <- normal_ID[grep("abnormal", normal_ID$SMPTHNTS),]
grep("Thyoird|Breast|Skin|Kidney|Lung", normal_ID_CD14$SMTS)




THCA_scatter$cancer <- 'THCA'
UCEC_scatter$cancer <- 'UCEC'
LUNG_scatter$cancer <- 'LUNG'
OV_scatter$cancer <- 'OV'
PRAD_scatter$cancer <- 'PRAD'
TNBC_yes_scatter$cancer <- 'TNBC_BRCA'
TNBC_no_scatter$cancer <- 'no_TNBC_BRCA'



scatter <- rbind(THCA_scatter,
                 UCEC_scatter,
                 LUNG_scatter,
                 OV_scatter,
                 PRAD_scatter,
                 TNBC_yes_scatter,
                 TNBC_no_scatter
                 )


ggscatter(scatter, x = "CD14", y = "ITGAM", size = 0.5, 
          rug = TRUE,                                
          color = "type", palette = "jco") +
  stat_cor(aes(color = type), method = "spearman")

library(plotly) #plotly 패키지 로드
plot_ly(scatter,x=~CD14,y=~ITGAM,z=~PDCD1,marker = list(size = 0.7, alpha = 0.5)) %>% add_markers(color=~type)


subset <- select(scatter, CD14,ITGAM,IL10,type)
colors <- c("#999999", "#E69F00", "#56B4E9")
colors <- colors[as.numeric(subset$type)]
s3d <- scatterplot3d(subset[,1:3], pch = 16, color=colors)
legend("bottom", legend = levels(scatter$type),
       col =  c("#999999", "#E69F00", "#56B4E9"), pch = 16, 
       inset = -0.25, xpd = TRUE, horiz = TRUE)



ggscatter(scatter, x = "CD14", y = c("SIRPA", "SIGLEC10","PDCD1"), size = 0.3,
          combine = TRUE, ylab = "Expression",
          color = "type", palette = "jco",
          add = "reg.line", conf.int = TRUE) +
  stat_cor(aes(color = type), method = "spearman", label.x = 4)



ggscatter(scatter, x = "CD14", y = "SIRPA", size = 0.3,
          color = "type", palette = "npg",
          facet.by = "type", #scales = "free_x",
          add = "reg.line", conf.int = TRUE) +
  stat_cor(aes(color = type), method = "spearman", label.y = 2)



library(cowplot) 
# Main plot
pmain <- ggplot(scatter, aes(x = CD14, y = SIRPA, color = type))+
  geom_point()+
  ggpubr::color_palette("jco") +
  stat_cor(aes(color = type), method = "spearman")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = scatter, aes(x = CD14, fill = type),
               alpha = 0.4, size = 0.2)+
  ggpubr::fill_palette("jco")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = scatter, aes(x = SIRPA, fill = type),
               alpha = 0.5, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("jco")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)


# Scatter plot colored by groups ("Species")
sp <- ggscatter(scatter, x = "CD14", y = "CD68",
                color = "type", palette = "jco",
                size = 3, alpha = 0.6)+
  border()  +
  stat_cor(aes(color = type), method = "spearman")                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(scatter, "CD14", fill = "type",
                   palette = "jco")
yplot <- ggdensity(scatter, "CD68", fill = "type", 
                   palette = "jco")+
  rotate()
# Cleaning the plots
sp <- sp + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
library(cowplot)
plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))

setwd('../../../')
test <- read.csv('TNBC_BRCA_merge.csv', header = T, sep = ',', row.names = 1)

test <- t(test)
test <- select(as.data.frame(test), TNBC.Subtype, CXCL16, IL8, AREG)
test$TNBC.Subtype <- ifelse(test$TNBC.Subtype=='no TNBC',0,1)
test$TNBC.Subtype <- as.factor(test$TNBC.Subtype)
fit <- lm(TNBC.Subtype ~ CXCL16 + I(CXCL16^2), data=test)
fit <- lm(AREG ~ TNBC.Subtype, data=test)
summary(fit)
