


library(ggplot2)
library(SingleCaseES)
library(WebPower)



# diversity responses
d_cloutier <- read.csv("Model-output/effects_diversity/Cloutier_diversity_RR.csv")
d_cloutier_amf <- read.csv("Model-output/effects_diversity/Cloutier_diversity_AMF_RR.csv")
d_cloutier_pathogen <- read.csv("Model-output/effects_diversity/Cloutier_diversity_pathogen_RR.csv")
d_cloutier_saprotroph <- read.csv("Model-output/effects_diversity/Cloutier_diversity_saprotroph_RR.csv")

d_song_fun <- read.csv("Model-output/effects_diversity/Song_fun-diversity_RR.csv")
d_song_amf <- read.csv("Model-output/effects_diversity/Song_fun_AMF-diversity_RR.csv")
d_song_pathogen <- read.csv("Model-output/effects_diversity/Song_fun_pathogen-diversity_RR.csv")
d_song_saprotroph <- read.csv("Model-output/effects_diversity/Song_fun_saprotroph-diversity_RR.csv")

d_strom <- read.csv("Model-output/effects_diversity/Strom_diversity_RR.csv")
d_strom_amf <- read.csv("Model-output/effects_diversity/Strom_AMF_diversity_RR.csv")
d_strom_pathogen <- read.csv("Model-output/effects_diversity/Strom_pathogen_diversity_RR.csv")
d_strom_pathogen2 <- d_strom_pathogen %>%
  dplyr::group_by(Comparison) %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(unique(.x))))
d_strom_pathogen2 <- as.data.frame(d_strom_pathogen2)
d_strom_pathogen2$X <- c(1,2)
d_strom_saprotroph <- read.csv("Model-output/effects_diversity/Strom_saprotroph_diversity_RR.csv")


# carbon responses
c_cloutier_OM <- read.csv("Model-output/effects_diversity/Cloutier_carbon-OM_RR.csv")
c_cloutier_POXC <- read.csv("Model-output/effects_diversity/Cloutier_carbon-POXC_RR.csv")

c_song <- read.csv("Model-output/effects_diversity/Song_carbon_RR.csv")

c_strom <- read.csv("Model-output/effects_diversity/Strom_carbon_RR.csv")


# diversity-carbon correlations
cor_cloutier <- read.csv("Model-output/effects_diversity/Cloutier_lmm-carbon.csv")
cor_cloutier_AMF <- read.csv("Model-output/effects_diversity/Cloutier_lmm-carbon_amf.csv")
cor_cloutier_pathogen <- read.csv("Model-output/effects_diversity/Cloutier_lmm-carbon_pathogen.csv")
cor_cloutier_saprotroph <- read.csv("Model-output/effects_diversity/Cloutier_lmm-carbon_saprotroph.csv")

cor_song_fun <- read.csv("Model-output/effects_diversity/Song_lmm-carbon_fun.csv")
cor_song_amf <- read.csv("Model-output/effects_diversity/Song_lmm-carbon_fun_amf.csv")
cor_song_pathogen <- read.csv("Model-output/effects_diversity/Song_lmm-carbon_fun_pathogen.csv")
cor_song_saprotroph <- read.csv("Model-output/effects_diversity/Song_lmm-carbon_fun_saprotroph.csv")

cor_strom <- read.csv("Model-output/effects_diversity/Strom_lmm-carbon.csv")
cor_strom_AMF <- read.csv("Model-output/effects_diversity/Strom_lmm-carbon_amf.csv")
cor_strom_pathogen <- read.csv("Model-output/effects_diversity/Strom_lmm-carbon_pathogen.csv")
cor_strom_saprotroph <- read.csv("Model-output/effects_diversity/Strom_lmm-carbon_saprotroph.csv")



# plot correlations
cors <- c(cor_cloutier$value_OM[which(cor_cloutier$statistic=="slope")],
          cor_cloutier_AMF$value_SOC[which(cor_cloutier_AMF$statistic=="slope")],
          cor_cloutier_pathogen$value_SOC[which(cor_cloutier_pathogen$statistic=="slope")],
          cor_cloutier_saprotroph$value_SOC[which(cor_cloutier_saprotroph$statistic=="slope")],

          cor_song_fun$value_SOM[which(cor_song_fun$statistic=="slope")],
          cor_song_amf$value_SOM[which(cor_song_amf$statistic=="slope")],
          cor_song_pathogen$value_SOM[which(cor_song_pathogen$statistic=="slope")],
          cor_song_saprotroph$value_SOM[which(cor_song_saprotroph$statistic=="slope")],
          
          cor_strom$value_SOC[which(cor_strom$statistic=="slope")],
          cor_strom_AMF$value_SOC[which(cor_strom_AMF$statistic=="slope")],
          cor_strom_pathogen$value_SOC[which(cor_strom_pathogen$statistic=="slope")],
          cor_strom_saprotroph$value_SOC[which(cor_strom_saprotroph$statistic=="slope")])


pvals <- c(cor_cloutier$value_OM[which(cor_cloutier$statistic=="p")],
          cor_cloutier_AMF$value_SOC[which(cor_cloutier_AMF$statistic=="p")],
          cor_cloutier_pathogen$value_SOC[which(cor_cloutier_pathogen$statistic=="p")],
          cor_cloutier_saprotroph$value_SOC[which(cor_cloutier_saprotroph$statistic=="p")],
          
          cor_song_fun$value_SOM[which(cor_song_fun$statistic=="p")],
          cor_song_amf$value_SOM[which(cor_song_amf$statistic=="p")],
          cor_song_pathogen$value_SOM[which(cor_song_pathogen$statistic=="p")],
          cor_song_saprotroph$value_SOM[which(cor_song_saprotroph$statistic=="p")],
          
          cor_strom$value_SOC[which(cor_strom$statistic=="p")],
          cor_strom_AMF$value_SOC[which(cor_strom_AMF$statistic=="p")],
          cor_strom_pathogen$value_SOC[which(cor_strom_pathogen$statistic=="p")],
          cor_strom_saprotroph$value_SOC[which(cor_strom_saprotroph$statistic=="p")])


names <- rep(c("Cloutier","Song","Strom"), each=4)
ns <- rep(c("N = 72","N = 12","N = 96"), each=4)

microbialgroup <- rep(c("All fungi", "AMF", "Pathogen", "Saprotroph"), 3)

cormat2 <- data.frame(names=names, 
                     microbialgroup=microbialgroup,
                     cors=cors,
                     pvals=pvals,
                     n=ns)
write.csv(cormat2, "Processed-data/effect-sizes/Betas_fun.csv")



library(plotrix)

# while plotrix is loaded anyway:
# set colors with color.scale
# need data as matrix*
cormat2a <- data.frame(corsAll=cormat2$cors[which(cormat2$microbialgroup=="All fungi")],
                       corsAMF=cormat2$cors[which(cormat2$microbialgroup=="AMF")],
                       corsPath=cormat2$cors[which(cormat2$microbialgroup=="Pathogen")],
                       corsSapr=cormat2$cors[which(cormat2$microbialgroup=="Saprotroph")])
cormat2b <- data.frame(corsAll=cormat2$pvals[which(cormat2$microbialgroup=="All fungi")],
                       corsAMF=cormat2$pvals[which(cormat2$microbialgroup=="AMF")],
                       corsPath=cormat2$pvals[which(cormat2$microbialgroup=="Pathogen")],
                       corsSapr=cormat2$pvals[which(cormat2$microbialgroup=="Saprotroph")])
cormat3 <- as.matrix(cormat2a)
rownames(cormat3) <- cormat2$names[c(1,5,9)]
colnames(cormat3) <- c("All", "AMF", "Pathogens", "Saprotrophs")

cormat3 <- cormat3[,c(1,2,4,3)]

colorfunc = colorRamp(c("red","white","blue"))
cormat3cols <- rgb(colorfunc((c(cormat3)+1)/2), maxColorValue = 255)

library(plotrix)
pdf("Figures/effect-plots/*3_diversity-carbon_betas_fungal-gropus.pdf", height=6, width=7)
par(mar = c(0.5, 8, 10, 2.1))
color2D.matplot(cormat3, 
                cellcolors=cormat3cols,
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black")
text(x = 1:4-0.6,
     y = par("usr")[3] + 3.1, adj=0,
     labels = colnames(cormat3),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 55,
     cex = 2)
axis(2, at = seq_len(nrow(cormat3)) -0.5,
     labels = rev(rownames(cormat3)), tick = FALSE, las = 1, cex.axis = 2)
axis(2, at = seq_len(nrow(cormat3)) -0.75,font=3,
     labels = rev(cormat2$n[c(1,5,9)]), tick = FALSE, las = 1, cex.axis = 1.6)
#points(x=1-0.15, y=3-0.4, pch=8)
#points(x=3-0.15, y=3-0.4, pch=8)
mtext("B. Fungi", side=3, line=6, cex=3,las=1, adj=0)
#mtext(expression(paste("A. Ashworth ",italic("et al.")," (2017) at MTREC")), side=3, adj=0, line=0.25, cex=1, font=4)

dev.off()



# adjust treatment names of cloutier
d_cloutier$Comparison <- gsub(" Spp", "-Spp", d_cloutier$Comparison)
c_cloutier_OM$Comparison <- gsub(" Spp", "-Spp", c_cloutier_OM$Comparison)


## set parameters for plotting
ptsize <- 5
ptshape <- 21
ptshapesig <- 19
ptbg <- "white"
boxthick <- 2
boxcol <- "black"
ylim <- rbind(c(-2,  2.6)) 
zerolinetype <- "solid"
zerolinecol <- "black"
zerolinethick <- 0.5
ebwid <- 0.1
addline_format <- function(x,...){gsub('\\s','\n',x)}
gridcol <- "lightgray"
txsize <- 6
stripsize <-12
axissize <- 11
axistext <- 10
loc_left <- margin(0.2, 0.1, 0, 0.2, "cm")
loc_middle <- margin(0.2, 0.1, 0, -0.1, "cm")
loc_right <- margin(0.2, 0.2, 0, 0, "cm")
second.axis <- sec_axis(~100* (exp(.)-1),name = "Percent change",                                          
                        breaks = c(-75, -50, 0, 50, 200,500,1000,2000),
                        labels = function(b) { paste0(round(b, 0))})
off_amf <- 0.15
off_sap <- 0.3
off_path <- 0.45
col_amf <- "cornflowerblue"
col_sap <- "darkgoldenrod4"
col_path <- "brown3"
# scale_y_continuous(sec.axis = second.axis)
#theme(axis.text.y.left = element_blank())
#theme(axis.text.y = element_blank())

# cloutier fung div
d_cloutier$Comparisonlab <- gsub(" ", "\n", d_cloutier$Comparison)
len <- dim(d_cloutier)[1]
d_cloutier$Comparison <- as.factor(d_cloutier$Comparison)
d_cloutier$Comparisonlab <- as.factor(d_cloutier$Comparisonlab)
d_cloutier$Comparison <- factor(d_cloutier$Comparison, levels(d_cloutier$Comparison)[c(3:8,1:2)])
d_cloutier$Comparisonlab <- factor(d_cloutier$Comparisonlab, levels(d_cloutier$Comparisonlab)[c(3:8,1:2)])
p_cloutier1 <- ggplot(dat=d_cloutier, aes(x=X, y=Est)) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal", plot.margin = loc_left, panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  labs(y="LRR", x="", title=expression(paste("A. Cloutier ",italic("et al.")," (2020)"))) +
  scale_x_continuous(breaks=c(1:8)+0.2,labels=d_cloutier$Comparisonlab, limits=c(1,8.5)) +
  scale_y_continuous(limits = ylim[1,], sec.axis = second.axis) +
# all fungi
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper),  width=0, size=1.5) +
  geom_point(data=d_cloutier, aes(x=X, y=Est, color="All"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_cloutier[which(sign(d_cloutier$CI_lower)==sign(d_cloutier$CI_upper)),], shape=ptshapesig, size=ptsize, aes(y=Est)) +
# amf
  geom_errorbar(data=d_cloutier_amf, color=col_amf, aes(ymin=CI_lower, ymax=CI_upper, x=X+off_amf), width=0, size=1.5) +
  geom_point(data=d_cloutier_amf, aes(x=X+off_amf, y=Est, col="AMF"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_cloutier_amf[which(sign(d_cloutier_amf$CI_lower)==sign(d_cloutier_amf$CI_upper)),], shape=ptshapesig, size=ptsize, col=col_amf, aes(y=Est, x=X+off_amf)) +
# sap
  geom_errorbar(data=d_cloutier_saprotroph, color=col_sap, aes(ymin=CI_lower, ymax=CI_upper, x=X+off_sap), width=0, size=1.5) +
  geom_point(data=d_cloutier_saprotroph, aes(x=X+off_sap, y=Est, col="Saprotrophs"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_cloutier_saprotroph[which(sign(d_cloutier_saprotroph$CI_lower)==sign(d_cloutier_saprotroph$CI_upper)),], shape=ptshapesig, size=ptsize, col=col_sap, aes(y=Est, x=X+off_sap)) +
# path
  geom_errorbar(data=d_cloutier_pathogen, color=col_path, aes(ymin=CI_lower, ymax=CI_upper, x=X+off_path), width=0, size=1.5) +
  geom_point(data=d_cloutier_pathogen, aes(x=X+off_path, y=Est, col="Pathogens"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_cloutier_pathogen[which(sign(d_cloutier_pathogen$CI_lower)==sign(d_cloutier_pathogen$CI_upper)),], shape=ptshapesig, size=ptsize, col=col_path, aes(y=Est, x=X+off_path)) +
  scale_color_manual(name = "Fungal group",
                     breaks = c("All", "AMF", "Saprotrophs", "Pathogens"),
                     values = c("All" = "black", "AMF" = col_amf, "Saprotrophs" = col_sap, "Pathogens" =  col_path) )
p_cloutier1




# song fung div
d_song_fun$Comparisonlab <- d_song_fun$Comparison
d_song_fun$Comparisonlab <- gsub("vs. ", "vs.-", d_song_fun$Comparisonlab)
d_song_fun$Comparisonlab <- gsub(" ", "\n", d_song_fun$Comparisonlab)
d_song_fun$Comparisonlab <- gsub("vs.-", "vs. ", d_song_fun$Comparisonlab)
len <- dim(d_song_fun)[1]
d_song_fun$Comparison <- as.factor(d_song_fun$Comparison)
d_song_fun$Comparisonlab <- as.factor(d_song_fun$Comparisonlab)
d_song_fun$Comparison <- factor(d_song_fun$Comparison, levels(d_song_fun$Comparison)[c(2,4,1,3)])
d_song_fun$Comparisonlab <- factor(d_song_fun$Comparisonlab, levels(d_song_fun$Comparisonlab)[c(2,4,1,3)])

p_song1 <- ggplot(dat=d_song_fun, aes(x=X, y=Est)) +
  theme_bw() +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=15), legend.position="bottom", legend.box = "horizontal", plot.margin = loc_left, panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  guides(colour = guide_legend(override.aes = list(size=5, shape=19))) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  labs(y="LRR", x="", title=expression(paste("B. Song ",italic("et al.")," (2018)"))) +
  scale_x_continuous(breaks=c(1:4)+0.2,labels=d_song_fun$Comparisonlab, limits=c(1,4.5)) +
  scale_y_continuous(limits = ylim[1,]) +
  # all fungi
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper),  width=0, size=1.5) +
  geom_point(data=d_song_fun, aes(x=X, y=Est, color="All"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_song_fun[which(sign(d_song_fun$CI_lower)==sign(d_song_fun$CI_upper)),], shape=ptshapesig, size=ptsize, aes(y=Est)) +
  # amf
  geom_errorbar(data=d_song_amf, color=col_amf, aes(ymin=CI_lower, ymax=CI_upper, x=X+off_amf), width=0, size=1.5) +
  geom_point(data=d_song_amf, aes(x=X+off_amf, y=Est, col="AMF"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_song_amf[which(sign(d_song_amf$CI_lower)==sign(d_song_amf$CI_upper)),], shape=ptshapesig, size=ptsize, col=col_amf, aes(y=Est, x=X+off_amf)) +
  # sap
  geom_errorbar(data=d_song_saprotroph, color=col_sap, aes(ymin=CI_lower, ymax=CI_upper, x=X+off_sap), width=0, size=1.5) +
  geom_point(data=d_song_saprotroph, aes(x=X+off_sap, y=Est, col="Saprotrophs"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_song_saprotroph[which(sign(d_song_saprotroph$CI_lower)==sign(d_song_saprotroph$CI_upper)),], shape=ptshapesig, size=ptsize, col=col_sap, aes(y=Est, x=X+off_sap)) +
  # path
  geom_errorbar(data=d_song_pathogen, color=col_path, aes(ymin=CI_lower, ymax=CI_upper, x=X+off_path), width=0, size=1.5) +
  geom_point(data=d_song_pathogen, aes(x=X+off_path, y=Est, col="Pathogens"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_song_pathogen[which(sign(d_song_pathogen$CI_lower)==sign(d_song_pathogen$CI_upper)),], shape=ptshapesig, size=ptsize, col=col_path, aes(y=Est, x=X+off_path)) +
  scale_color_manual(name = "Fungal group",
                     breaks = c("All", "AMF", "Saprotrophs", "Pathogens"),
                     values = c("All" = "black", "AMF" = col_amf, "Saprotrophs" = col_sap, "Pathogens" =  col_path) )
p_song1





# strom fun diversity
d_strom$Comparisonlab <- gsub(" ", "\n", d_strom$Comparison)
len <- dim(d_strom)[1]
d_strom$Comparisonlab <- as.factor(d_strom$Comparisonlab)
p_strom1 <- ggplot(dat=d_strom, aes(X, y=Est)) +
  theme_bw() +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=15), legend.position="bottom", legend.box = "horizontal", plot.margin = loc_left, axis.text.y.left = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  guides(colour = guide_legend(override.aes = list(size=5, shape=19))) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  labs(y="", x="", title=expression(paste("C. Strom ",italic("et al.")," (2020)"))) +
  scale_x_continuous(breaks=c(1,2)+0.2,labels=d_strom$Comparisonlab, limits=c(0.75,2.75)) +
  scale_y_continuous(limits = ylim[1,], sec.axis = second.axis) +
  # all fungi
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper),  width=0, size=1.5) +
  geom_point(data=d_strom, aes(x=X, y=Est, color="All"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_strom[which(sign(d_strom$CI_lower)==sign(d_strom$CI_upper)),], shape=ptshapesig, size=ptsize, aes(y=Est)) +
  # amf
  geom_errorbar(data=d_strom_amf, color=col_amf, aes(ymin=CI_lower, ymax=CI_upper, x=X+off_amf), width=0, size=1.5) +
  geom_point(data=d_strom_amf, aes(x=X+off_amf, y=Est, col="AMF"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_strom_amf[which(sign(d_strom_amf$CI_lower)==sign(d_strom_amf$CI_upper)),], shape=ptshapesig, size=ptsize, col=col_amf, aes(y=Est, x=X+off_amf)) +
  # sap
  geom_errorbar(data=d_strom_saprotroph, color=col_sap, aes(ymin=CI_lower, ymax=CI_upper, x=X+off_sap), width=0, size=1.5) +
  geom_point(data=d_strom_saprotroph, aes(x=X+off_sap, y=Est, col="Saprotrophs"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_strom_saprotroph[which(sign(d_strom_saprotroph$CI_lower)==sign(d_strom_saprotroph$CI_upper)),], shape=ptshapesig, size=ptsize, col=col_sap, aes(y=Est, x=X+off_sap)) +
  # path
  geom_errorbar(data=d_strom_pathogen2, color=col_path, aes(ymin=CI_lower, ymax=CI_upper, x=X+off_path), width=0, size=1.5) +
  geom_point(data=d_strom_pathogen2, aes(x=X+off_path, y=Est, col="Pathogens"), size=ptsize, shape=ptshape, fill=ptbg) +
  geom_point(data=d_strom_pathogen2[which(sign(d_strom_pathogen2$CI_lower)==sign(d_strom_pathogen2$CI_upper)),], shape=ptshapesig, size=ptsize, col=col_path, aes(y=Est, x=X+off_path)) +
  scale_color_manual(name = "Fungal group",
                     breaks = c("All", "AMF", "Saprotrophs", "Pathogens"),
                     values = c("All" = "black", "AMF" = col_amf, "Saprotrophs" = col_sap, "Pathogens" =  col_path) )
p_strom1







p_song_strom <- ggpubr::ggarrange(p_song1,  p_strom1, ncol=2, widths = c(1,0.7), common.legend = TRUE, legend="bottom")
p_css <- ggpubr::ggarrange(p_cloutier1,  p_song_strom, ncol=1, nrow=2)

pdf("Figures/effect-plots/*1_LRR_fungal-functional-groups.pdf", height=7, width=10)
p_css
dev.off()




