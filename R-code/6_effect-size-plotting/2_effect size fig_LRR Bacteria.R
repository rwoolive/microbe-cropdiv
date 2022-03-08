


library(ggplot2)
library(SingleCaseES)
library(WebPower)
library(dplyr)



# diversity responses
d_ashworth_mtrec_cover <- read.csv("Model-output/effects_diversity/Ashworth_diversity-mtrec-cover_RR.csv")
d_ashworth_mtrec_cropsys <- read.csv("Model-output/effects_diversity/Ashworth_diversity-mtrec-cropsys_RR.csv")
d_ashworth_recm_cover <- read.csv("Model-output/effects_diversity/Ashworth_diversity-recm-cover-yr_RR.csv")
d_ashworth_recm_cropsys <- read.csv("Model-output/effects_diversity/Ashworth_diversity-recm-cropsys_RR.csv")
d_gao <- read.csv("Model-output/effects_diversity/Gao_diversity_RR.csv")
d_song_bac <- read.csv("Model-output/effects_diversity/Song_bac-diversity_RR.csv")


# # carbon responses
# c_ashworth_mtrec_cover <- read.csv("Model-output/effects_diversity/Ashworth_carbon-mtrec-cover_RR.csv")
# c_ashworth_mtrec_cropsys <- read.csv("Model-output/effects_diversity/Ashworth_carbon-mtrec-cropsys_RR.csv")
# c_ashworth_recm_cover <- read.csv("Model-output/effects_diversity/Ashworth_carbon-recm-cover_RR.csv")
# c_ashworth_recm_cropsys <- read.csv("Model-output/effects_diversity/Ashworth_carbon-recm-cropsys_RR.csv")
# c_gao_DOC <- read.csv("Model-output/effects_diversity/Gao_carbon-DOC_RR.csv")
# c_gao_SOC <- read.csv("Model-output/effects_diversity/Gao_carbon-SOC_RR.csv")
# c_song <- read.csv("Model-output/effects_diversity/Song_carbon_RR.csv")









# diversity-carbon correlations
cor_ashworth <- read.csv("Model-output/effects_diversity/Ashworth_lmm-carbon.csv")
cor_gao <- read.csv("Model-output/effects_diversity/Gao_lmm-carbon.csv")
cor_song <- read.csv("Model-output/effects_diversity/Song_lmm-carbon_bac.csv")



# plot correlations
cors <- c(cor_ashworth$value_mtrec[which(cor_ashworth$statistic=="slope")],
          cor_ashworth$value_recm[which(cor_ashworth$statistic=="slope")],
          cor_gao$value_SOC[which(cor_gao$statistic=="slope")],
          cor_song$value_SOM[which(cor_song$statistic=="slope")])

pvals <- c(cor_ashworth$value_mtrec[which(cor_ashworth$statistic=="p")],
           cor_ashworth$value_recm[which(cor_ashworth$statistic=="p")],
           cor_gao$value_SOC[which(cor_gao$statistic=="p")],
           cor_song$value_SOM[which(cor_song$statistic=="p")])

names <- rep(c("Ashworth\nat MTREC","Ashworth\nat RECM","Gao","Song"), each=1)
ns <- rep(c("N = 90","N = 96","N = 12","N = 12"), each=1)

microbialgroup <- rep(c("All"), 4)

cormat2 <- data.frame(names=names, 
                      microbialgroup=microbialgroup,
                      cors=cors,
                      pvals=pvals,
                      n=ns)
write.csv(cormat2, "Processed-data/effect-sizes/Betas-bac.csv")



library(plotrix)

# while plotrix is loaded anyway:
# set colors with color.scale
# need data as matrix*
cormat2a <- data.frame(corsAll=cormat2$cors[which(cormat2$microbialgroup=="All")])
cormat2b <- data.frame(corsAll=cormat2$pvals[which(cormat2$microbialgroup=="All")])
cormat3 <- as.matrix(cormat2a)
rownames(cormat3) <- cormat2$names[c(1:4)]
colnames(cormat3) <- c("All")


colorfunc = colorRamp(c("red","white","blue"))
cormat3cols <- rgb(colorfunc((c(cormat3)+1)/2), maxColorValue = 255)


library(plotrix)
pdf("Figures/effect-plots/*3_diversity-carbon_betas_bacteria.pdf", height=7, width=4)
par(mar = c(0.5, 11, 10, 2.1))
color2D.matplot(cormat3, 
                cellcolors = cormat3cols,
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black")
text(x = 1:4-0.6,
     y = par("usr")[3] + 4.2, adj=0,
     labels = colnames(cormat3),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 55,
     cex = 2)
axis(2, at = seq_len(nrow(cormat3))-c(0.5,0.5,0.35,0.35),
     labels = rev(rownames(cormat3)), tick = FALSE, las = 1, cex.axis = 2)
axis(2, at = seq_len(nrow(cormat3)) -0.75,font=3,
     labels = rev(cormat2$n), tick = FALSE, las = 1, cex.axis = 1.6)
#text(x=1-0.2, y=1-0.6, ".", cex=3)
mtext("A. Bacteria", sisde=3, line=6, cex=3,las=1, at = 0)#adj=0
#mtext(expression(paste("A. Ashworth ",italic("et al.")," (2017) at MTREC")), side=3, adj=0, line=0.25, cex=1, font=4)
dev.off()

colorfunc = colorRamp(c("blue","white","red"))
vals <- as.matrix(seq(-1,1,length.out = 100)+1)/2
allcols <- rgb(colorfunc(vals), maxColorValue = 255)
pdf("Figures/effect-plots/*3_diversity-carbon_betas_scale.pdf", height=1, width=3)
par(mar = c(1, 1, 1, 1))
legend_image <- as.raster(matrix(allcols))
plot(c(-1,1),c(-1.5,2),type = 'n', axes = F,xlab = '', ylab = '', main = 'Beta')
text(y=-1, x = seq(-1,1,l=5), labels =round( seq(-1,1,l=5),1))
rasterImage(t(legend_image), -1, -0.5, 1,1.5)
dev.off()





## Not run: 
# requires viridisLite
library(viridisLite)
plot(0,xlim=c(-1,1),ylim=c(-1,1),type="n",axes=FALSE,
     main="Approximating other color scales",xlab="",ylab="")
gradient.rect(-1,0.8,1,0.95,nslices=50,
              col=color.scale(1:50,1,
                              rev(c(0,0.3,0.6,0.8,1,1)),
                             rev( c(0,0,0,0,0,0,1))))
text(0,1,"color.scale")
gradient.rect(-1,0.65,1,0.8,col=heat.colors(50))
text(0,0.6,"heat.colors")
gradient.rect(-1,0.3,1,0.45,nslices=50,
              col=color.scale(1:50,c(0,0.2,0.9,0.95,0.95),
                              c(0.7,0.8,0.9,0.7,0.95),
                              c(0.1,0,0,0.35,0.95)))
text(0,0.5,"color.scale")
gradient.rect(-1,0.15,1,0.3,col=terrain.colors(50))
text(0,0.1,"terrain.colors")
gradient.rect(-1,-0.2,1,-0.05,nslices=50,
              col=color.scale(1:50,c(0.3,0,0.3,0.1,1,0.95,1),
                              c(0,0.3,0.9,1,1,0.85,0.85),
                              c(1,1,0.9,0.1,0,0.5,0.5)))
text(0,0,"color.scale")
gradient.rect(-1,-0.35,1,-0.2,col=topo.colors(50))
text(0,-0.4,"topo.colors")
gradient.rect(-1,-0.7,1,-0.55,nslices=50,
              col=color.scale(1:50,c(0.3,0.2,0,0.4,0.95),
                              c(0.1,0.3,0.6,0.75,0.95),
                              c(0.3,0.6,0.5,0.4,0)))
text(0,-0.5,"color.scale")
gradient.rect(-1,-0.85,1,-0.7,col=viridis(50))
text(0,-0.9,"viridis")







## set parameters for plotting
ptsize <- 3
ptshape <- 21
ptshapesig <- 19
ptbg <- "white"
boxthick <- 2
boxcol <- "black"
ylim <- rbind(c(-1,  4.25)) 
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
loc_whole <- margin(0.2, 0.2, 0, 0.2, "cm")
loc_right <- margin(0.2, 0.2, 0, 0, "cm")
second.axis <- sec_axis(~100* (exp(.)-1),name = "Percent change",                                          
                        breaks = c( -50, 0, 50,500,2000,5000),
                        labels = function(b) { paste0(round(b, 0))})


# mtrec bac diversity by cover
d_ashworth_mtrec_cover$Comparisonlab <- gsub(" ", "\n", d_ashworth_mtrec_cover$Comparison)
len <- dim(d_ashworth_mtrec_cover)[1]
p_ashworth1 <- ggplot(dat=d_ashworth_mtrec_cover, aes(x=Comparison, y=Est)) +
  theme_bw() +
  theme(plot.margin = loc_left, panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="LRR", x="", title="A") +
  scale_x_discrete(labels =d_ashworth_mtrec_cover$Comparisonlab) +
  geom_point(data=d_ashworth_mtrec_cover[which(sign(d_ashworth_mtrec_cover$CI_lower)==sign(d_ashworth_mtrec_cover$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))+
  scale_y_continuous(limits = ylim[1,]) 
p_ashworth1

# mtrec bac diversity by cropsys
d_ashworth_mtrec_cropsys$Comparisonlab <- gsub(" ", "\n", d_ashworth_mtrec_cropsys$Comparison)
len <- dim(d_ashworth_mtrec_cropsys)[1]
p_ashworth2 <- ggplot(dat=d_ashworth_mtrec_cropsys, aes(x=Comparison, y=Est)) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal", plot.margin = loc_right, axis.text.y.left = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="", x="", title="B") +
  scale_x_discrete(labels =d_ashworth_mtrec_cropsys$Comparisonlab)+
  geom_point(data=d_ashworth_mtrec_cropsys[which(sign(d_ashworth_mtrec_cropsys$CI_lower)==sign(d_ashworth_mtrec_cropsys$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est)) +
  scale_y_continuous(limits = ylim[1,], sec.axis = second.axis) 
p_ashworth2


# # mtrec carbon by cover
# c_ashworth_mtrec_cover$Comparisonlab <- gsub(" ", "\n", c_ashworth_mtrec_cover$Comparison)
# len <- dim(c_ashworth_mtrec_cover)[1]
# p_ashworth3 <- ggplot(dat=c_ashworth_mtrec_cover, aes(x=Comparison, y=Est)) +
#   theme_bw() +
#   theme(plot.margin = margin(0.2, 0, 0, -0.1, "cm"), axis.text.y = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
#   geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
#   geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
#   geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
#   labs(y="", x="", title="C") +
#   scale_x_discrete(labels =c_ashworth_mtrec_cover$Comparisonlab)+
#   geom_point(data=c_ashworth_mtrec_cover[which(sign(c_ashworth_mtrec_cover$CI_lower)==sign(c_ashworth_mtrec_cover$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))+
#   scale_y_continuous(limits = ylim[1,])
# p_ashworth3
# 
# # mtrec carbon by cropsys
# c_ashworth_mtrec_cropsys$Comparisonlab <- gsub(" ", "\n", c_ashworth_mtrec_cropsys$Comparison)
# len <- dim(c_ashworth_mtrec_cropsys)[1]
# p_ashworth4 <- ggplot(dat=c_ashworth_mtrec_cropsys, aes(x=Comparison, y=Est)) +
#   theme_bw() +
#   theme(plot.margin = margin(0.2, 0.2, 0, 0, "cm"), axis.text.y.left = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
#   geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
#   geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
#   geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
#   labs(y="", x="", title="D") +
#   scale_x_discrete(labels =c_ashworth_mtrec_cropsys$Comparisonlab)+
#   geom_point(data=c_ashworth_mtrec_cropsys[which(sign(c_ashworth_mtrec_cropsys$CI_lower)==sign(c_ashworth_mtrec_cropsys$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est)) +
#   scale_y_continuous(limits = ylim[1,], sec.axis = sec_axis(~100* (exp(.)-1),name = "Percent change",                                          
#                                                             breaks = c(-50, 0, 50, 100, 200,500,1000,2000),
#                                                             labels = function(b) { paste0(round(b, 0))})) 
# p_ashworth4


# recm bac diversity by cover (this one has interaction with year)
d_ashworth_recm_cover$Year <- c(2013, 2014, 2013, 2014)
d_ashworth_recm_cover$Comparisonlab <- gsub(" ", "\n", d_ashworth_recm_cover$Comparison)
len <- dim(d_ashworth_recm_cover)[1]
p_ashworth5 <- ggplot(dat=d_ashworth_recm_cover, aes(x=Comparison, y=Est)) +
  theme_bw() +
  theme(plot.margin = loc_left, panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="LRR", x="", title="C") +
  scale_x_discrete(labels =d_ashworth_recm_cover$Comparisonlab[c(1,3)]) +
  facet_grid(.~Year, labeller = label_both) +
  geom_point(data=d_ashworth_recm_cover[which(sign(d_ashworth_recm_cover$CI_lower)==sign(d_ashworth_recm_cover$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))+
  scale_y_continuous(limits = ylim[1,]) 
p_ashworth5

# recm bac diversity by cropsys
d_ashworth_recm_cropsys$Comparisonlab <- gsub(" ", "\n", d_ashworth_recm_cropsys$Comparison)
len <- dim(d_ashworth_recm_cropsys)[1]
p_ashworth6 <- ggplot(dat=d_ashworth_recm_cropsys, aes(x=Comparison, y=Est)) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal", plot.margin = loc_right, axis.text.y.left = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="", x="", title="D") +
  scale_x_discrete(labels =d_ashworth_recm_cropsys$Comparisonlab)+
  geom_point(data=d_ashworth_recm_cropsys[which(sign(d_ashworth_recm_cropsys$CI_lower)==sign(d_ashworth_recm_cropsys$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))+
  scale_y_continuous(limits = ylim[1,], sec.axis = second.axis) 
p_ashworth6




# # recm carbon by cover
# c_ashworth_recm_cover$Comparisonlab <- gsub(" ", "\n", c_ashworth_recm_cover$Comparison)
# len <- dim(c_ashworth_recm_cover)[1]
# p_ashworth7 <- ggplot(dat=c_ashworth_recm_cover, aes(x=Comparison, y=Est)) +
#   theme_bw() +
#   theme(plot.margin = margin(0.2, 0, 0, -0.1, "cm"), axis.text.y = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
#   geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
#   geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
#   geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
#   labs(y="", x="", title="G") +
#   scale_x_discrete(labels =c_ashworth_recm_cover$Comparisonlab)+
#   geom_point(data=c_ashworth_recm_cover[which(sign(c_ashworth_recm_cover$CI_lower)==sign(c_ashworth_recm_cover$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))+
#   scale_y_continuous(limits = ylim[1,])
# p_ashworth7
# 
# # recm carbon by cropsys
# c_ashworth_recm_cropsys$Comparisonlab <- gsub(" ", "\n", c_ashworth_recm_cropsys$Comparison)
# len <- dim(c_ashworth_recm_cropsys)[1]
# p_ashworth8 <- ggplot(dat=c_ashworth_recm_cropsys, aes(x=Comparison, y=Est)) +
#   theme_bw() +
#   theme(plot.margin = margin(0.2, 0.2, 0, 0, "cm"), axis.text.y.left = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
#   geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
#   geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
#   geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
#   labs(y="", x="", title="H") +
#   scale_x_discrete(labels =c_ashworth_recm_cropsys$Comparisonlab)+
#   geom_point(data=c_ashworth_recm_cropsys[which(sign(c_ashworth_recm_cropsys$CI_lower)==sign(c_ashworth_recm_cropsys$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est)) +
#   scale_y_continuous(limits = ylim[1,], sec.axis = sec_axis(~100* (exp(.)-1),name = "Percent change",                                          
#                                                             breaks = c(-50, 0, 50, 100, 200,500,1000,2000),
#                                                             labels = function(b) { paste0(round(b, 0))})) 
# p_ashworth8


# gao bac
d_gao$Comparisonlab <- gsub(" ", "\n", d_gao$Comparison)
len <- dim(d_gao)[1]
d_gao$Comparison <- as.factor(d_gao$Comparison)
d_gao$Comparisonlab <- as.factor(d_gao$Comparisonlab)
d_gao$Comparison <- factor(d_gao$Comparison, levels(d_gao$Comparison)[c(2,3,1)])
d_gao$Comparisonlab <- factor(d_gao$Comparisonlab, levels(d_gao$Comparisonlab)[c(2,3,1)])
p_gao1 <- ggplot(dat=d_gao, aes(x=Comparison, y=Est)) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal", plot.margin = loc_whole, panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="LRR", x="", title="E") +
  scale_x_discrete(labels =d_gao$Comparisonlab)+
  geom_point(data=d_gao[which(sign(d_gao$CI_lower)==sign(d_gao$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))+
  scale_y_continuous(limits = ylim[1,], sec.axis = second.axis)
p_gao1



# song bac
d_song_bac$Comparisonlab <- gsub(" ", "\n", d_song_bac$Comparison)
len <- dim(d_song_bac)[1]
d_song_bac$Comparison <- as.factor(d_song_bac$Comparison)
d_song_bac$Comparisonlab <- as.factor(d_song_bac$Comparisonlab)
d_song_bac$Comparison <- factor(d_song_bac$Comparison, levels(d_song_bac$Comparison)[c(2,4,1,3)])
d_song_bac$Comparisonlab <- factor(d_song_bac$Comparisonlab, levels(d_song_bac$Comparisonlab)[c(2,4,1,3)])
p_song1 <- ggplot(dat=d_song_bac, aes(x=Comparison, y=Est)) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal", plot.margin = loc_whole, panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="LRR", x="", title="F") +
  scale_x_discrete(labels =d_song_bac$Comparisonlab)+
  geom_point(data=d_song_bac[which(sign(d_song_bac$CI_lower)==sign(d_song_bac$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))+
  scale_y_continuous(limits = ylim[1,], sec.axis = second.axis)
p_song1




# 
# 
# 
# # column titles
# d=data.frame(x1=c(-1.9), x2=c(1.9), y1=c(-10), y2=c(10))
# bacdiv <- ggplot() + lims(y=c(-10,10), x=c(-2,2)) +
#   geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="gray", color="black", alpha=0.5) + 
#   annotate(geom="text", x=0, y=0, label="Bacterial diversity", size=txsize) + theme_void() 
# carb <- ggplot() + lims(y=c(-10,10), x=c(-2,2)) +
#   geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="gray", color="black", alpha=0.5) +
#   annotate(geom="text", x=0, y=0, label="Carbon", size=txsize) + theme_void() 



# row titles
tash_m <- ggplot() + lims(y=c(-10,10), x=c(-2,2)) +
  annotate(geom="text", x=-0.75, y=0, label=expression(paste("Ashworth ",italic("et al.")," (2017)")), size=txsize, angle = 90) + 
  annotate(geom="text", x=0.75, y=0, label=expression(paste("MTREC")), size=txsize, angle = 90) +
  theme_void()
tash_r <- ggplot() + lims(y=c(-10,10), x=c(-2,2)) +
  annotate(geom="text", x=-0.75, y=0, label=expression(paste("Ashworth ",italic("et al.")," (2017)")), size=txsize, angle = 90) + 
  annotate(geom="text", x=0.75, y=0, label=expression(paste("RECM")), size=txsize, angle = 90) +
  theme_void()
tgao <- ggplot() + lims(y=c(-10,10), x=c(-2,2)) +
  annotate(geom="text", x=0, y=0, label=expression(paste("Gao ",italic("et al.")," (2019)")), size=txsize, angle = 90) + 
  theme_void()
tsong <- ggplot() + lims(y=c(-10,10), x=c(-2,2)) +
  annotate(geom="text", x=0, y=0, label=expression(paste("Song ",italic("et al.")," (2018)")), size=txsize, angle = 90) + 
  theme_void()

# NA
plot_NA <- ggplot() + lims(y=c(-10,10), x=c(-2,2)) +
  annotate(geom="text", x=0, y=0, label="NA", size=5) + theme_void()

# plot nothing
plot_no <- ggplot() + lims(y=c(-10,10), x=c(-2,2)) + theme_void()






pdf("Figures/effect-plots/*2_LRR-bacteria.pdf", height=12, width=8)
laymat <- rbind(c(1,2,2,3,3),
            c(4,5,5,6,6),
            c(7,8,8,8,8),
            c(9,10,10,10,10))
gridExtra::grid.arrange(tash_m, p_ashworth1, p_ashworth2, 
                        tash_r, p_ashworth5, p_ashworth6, 
                        tgao, p_gao1,
                        tsong, p_song1,
                  nrow=4,
                  ncol=5, widths=c(0.55,1,1,1,1),
                  layout_matrix=laymat)
dev.off()








