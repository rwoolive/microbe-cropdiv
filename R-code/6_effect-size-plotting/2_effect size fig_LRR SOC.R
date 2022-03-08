


library(ggplot2)
library(SingleCaseES)
library(WebPower)




# carbon responses
c_ashworth_mtrec_cover <- read.csv("Model-output/effects_diversity/Ashworth_carbon-mtrec-cover_RR.csv")
c_ashworth_mtrec_cropsys <- read.csv("Model-output/effects_diversity/Ashworth_carbon-mtrec-cropsys_RR.csv")
c_ashworth_recm_cover <- read.csv("Model-output/effects_diversity/Ashworth_carbon-recm-cover_RR.csv")
c_ashworth_recm_cropsys <- read.csv("Model-output/effects_diversity/Ashworth_carbon-recm-cropsys_RR.csv")
c_cloutier_OM <- read.csv("Model-output/effects_diversity/Cloutier_carbon-OM_RR.csv")
c_cloutier_POXC <- read.csv("Model-output/effects_diversity/Cloutier_carbon-POXC_RR.csv")
c_gao_DOC <- read.csv("Model-output/effects_diversity/Gao_carbon-DOC_RR.csv")
c_gao_SOC <- read.csv("Model-output/effects_diversity/Gao_carbon-SOC_RR.csv")
c_song <- read.csv("Model-output/effects_diversity/Song_carbon_RR.csv")
c_strom <- read.csv("Model-output/effects_diversity/Strom_carbon_RR.csv")



# adjust treatment names of cloutier
c_cloutier_OM$Comparison <- gsub(" Spp", "-Spp", c_cloutier_OM$Comparison)


## set parameters for plotting
## set parameters for plotting
ptsize <- 5
ptshape <- 21
ptshapesig <- 19
ptbg <- "white"
boxthick <- 2
boxcol <- "black"
ylim <- rbind(c(-0.52,  0.5)) 
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
                        breaks = c(-40, -25, 0, 10, 25, 50),
                        labels = function(b) { paste0(round(b, 0))})



# mtrec carbon by cover
c_ashworth_mtrec_cover$Comparisonlab <- gsub(" ", "\n", c_ashworth_mtrec_cover$Comparison)
len <- dim(c_ashworth_mtrec_cover)[1]
p_ashworth3 <- ggplot(dat=c_ashworth_mtrec_cover, aes(x=Comparison, y=Est)) +
  scale_y_continuous(limits = ylim[1,]) +
  theme_bw() +
  theme(plot.margin = loc_left, panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="LRR", x="", title="A") +
  scale_x_discrete(labels =c_ashworth_mtrec_cover$Comparisonlab)+
  geom_point(data=c_ashworth_mtrec_cover[which(sign(c_ashworth_mtrec_cover$CI_lower)==sign(c_ashworth_mtrec_cover$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))
p_ashworth3

# mtrec carbon by cropsys
c_ashworth_mtrec_cropsys$Comparisonlab <- gsub(" ", "\n", c_ashworth_mtrec_cropsys$Comparison)
len <- dim(c_ashworth_mtrec_cropsys)[1]
p_ashworth4 <- ggplot(dat=c_ashworth_mtrec_cropsys, aes(x=Comparison, y=Est)) +
  scale_y_continuous(limits = ylim[1,]) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal", plot.margin = loc_middle, axis.text.y = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="", x="", title="B") +
  scale_x_discrete(labels =c_ashworth_mtrec_cropsys$Comparisonlab)+
  geom_point(data=c_ashworth_mtrec_cropsys[which(sign(c_ashworth_mtrec_cropsys$CI_lower)==sign(c_ashworth_mtrec_cropsys$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))
p_ashworth4



# recm carbon by cover
c_ashworth_recm_cover$Comparisonlab <- gsub(" ", "\n", c_ashworth_recm_cover$Comparison)
len <- dim(c_ashworth_recm_cover)[1]
p_ashworth7 <- ggplot(dat=c_ashworth_recm_cover, aes(x=Comparison, y=Est)) +
  scale_y_continuous(limits = ylim[1,]) +
  theme_bw() +
  theme(plot.margin = loc_left, axis.text.y = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="", x="", title="C") +
  scale_x_discrete(labels =c_ashworth_recm_cover$Comparisonlab)+
  geom_point(data=c_ashworth_recm_cover[which(sign(c_ashworth_recm_cover$CI_lower)==sign(c_ashworth_recm_cover$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))
p_ashworth7

# recm carbon by cropsys
c_ashworth_recm_cropsys$Comparisonlab <- gsub(" ", "\n", c_ashworth_recm_cropsys$Comparison)
len <- dim(c_ashworth_recm_cropsys)[1]
p_ashworth8 <- ggplot(dat=c_ashworth_recm_cropsys, aes(x=Comparison, y=Est)) +
  scale_y_continuous(limits = ylim[1,], sec.axis = second.axis) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal", plot.margin = loc_right,  axis.text.y.left = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="", x="", title="D") +
  scale_x_discrete(labels =c_ashworth_recm_cropsys$Comparisonlab)+
  geom_point(data=c_ashworth_recm_cropsys[which(sign(c_ashworth_recm_cropsys$CI_lower)==sign(c_ashworth_recm_cropsys$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))
p_ashworth8







# cloutier som
c_cloutier_OM$Comparisonlab <- gsub(" ", "\n", c_cloutier_OM$Comparison)
c_cloutier_OM$Comparisonlab <- gsub("mix\nvs.", "mix vs.", c_cloutier_OM$Comparisonlab)
len <- dim(c_cloutier_OM)[1]
c_cloutier_OM$Comparison <- as.factor(c_cloutier_OM$Comparison)
c_cloutier_OM$Comparisonlab <- as.factor(c_cloutier_OM$Comparisonlab)
c_cloutier_OM$Comparison <- factor(c_cloutier_OM$Comparison, levels(c_cloutier_OM$Comparison)[c(3:8,1:2)])
c_cloutier_OM$Comparisonlab <- factor(c_cloutier_OM$Comparisonlab, levels(c_cloutier_OM$Comparisonlab)[c(3:8,1:2)])
p_cloutier2 <- ggplot(dat=c_cloutier_OM, aes(x=Comparison, y=Est)) +
  scale_y_continuous(limits = ylim[1,]) +
  theme_bw() +
  theme(plot.margin = loc_left, panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="LRR", x="", title="E") +
  scale_x_discrete(labels =c_cloutier_OM$Comparisonlab)+
  geom_point(data=c_cloutier_OM[which(sign(c_cloutier_OM$CI_lower)==sign(c_cloutier_OM$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))
p_cloutier2


# cloutier poxc
c_cloutier_POXC$Comparisonlab <- gsub(" ", "\n", c_cloutier_POXC$Comparison)
len <- dim(c_cloutier_POXC)[1]
c_cloutier_POXC$Comparison <- as.factor(c_cloutier_POXC$Comparison)
c_cloutier_POXC$Comparisonlab <- as.factor(c_cloutier_POXC$Comparisonlab)
c_cloutier_POXC$Comparison <- factor(c_cloutier_POXC$Comparison, levels(c_cloutier_POXC$Comparison)[c(3:8,1:2)])
c_cloutier_POXC$Comparisonlab <- factor(c_cloutier_POXC$Comparisonlab, levels(c_cloutier_POXC$Comparisonlab)[c(3:8,1:2)])
p_cloutier3 <- ggplot(dat=c_cloutier_POXC, aes(x=Comparison, y=Est)) +
  scale_y_continuous(limits = ylim[1,]) +
  theme_bw() +
  theme(plot.margin = loc_left, panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="LRR", x="", title="E2") +
  scale_x_discrete(labels =c_cloutier_POXC$Comparisonlab)+
  geom_point(data=c_cloutier_POXC[which(sign(c_cloutier_POXC$CI_lower)==sign(c_cloutier_POXC$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est))
p_cloutier3






# gao doc
c_gao_DOC$Comparisonlab <- gsub(" ", "\n", c_gao_DOC$Comparison)
len <- dim(c_gao_DOC)[1]
c_gao_DOC$Comparison <- as.factor(c_gao_DOC$Comparison)
c_gao_DOC$Comparisonlab <- as.factor(c_gao_DOC$Comparisonlab)
c_gao_DOC$Comparison <- factor(c_gao_DOC$Comparison, levels(c_gao_DOC$Comparison)[c(2,3,1)])
c_gao_DOC$Comparisonlab <- factor(c_gao_DOC$Comparisonlab, levels(c_gao_DOC$Comparisonlab)[c(2,3,1)])
p_gao2 <- ggplot(dat=c_gao_DOC, aes(x=Comparison, y=Est)) +
  scale_y_continuous(limits = ylim[1,], sec.axis = second.axis) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal", plot.margin = loc_right,  axis.text.y.left = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  #geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="", x="", title="F2") +
  scale_x_discrete(labels =c_gao_DOC$Comparisonlab) +
  geom_point(data=c_gao_DOC[1:3,], shape=ptshapesig, size=ptsize, aes(x=X, y=Est)) # intercropping reduced DOC for all treatments (see Gao 2019)
p_gao2



# gao soc 
c_gao_SOC$Comparisonlab <- gsub(" ", "\n", c_gao_SOC$Comparison)
len <- dim(c_gao_SOC)[1]
c_gao_SOC$Comparison <- as.factor(c_gao_SOC$Comparison)
c_gao_SOC$Comparisonlab <- as.factor(c_gao_SOC$Comparisonlab)
c_gao_SOC$Comparison <- factor(c_gao_SOC$Comparison, levels(c_gao_SOC$Comparison)[c(2,3,1)])
c_gao_SOC$Comparisonlab <- factor(c_gao_SOC$Comparisonlab, levels(c_gao_SOC$Comparisonlab)[c(2,3,1)])
p_gao3 <- ggplot(dat=c_gao_SOC, aes(x=Comparison, y=Est)) +
  scale_y_continuous(limits = ylim[1,], sec.axis = second.axis) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal", plot.margin = loc_right,  axis.text.y.left = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  #geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="", x="", title="F") +
  scale_x_discrete(labels =c_gao_SOC$Comparisonlab) # no significant effect of age on SOC (see Gao 2019)

p_gao3





# song carbon
c_song$Comparisonlab <- gsub(" ", "\n", c_song$Comparison)
c_song$Comparisonlab <- gsub("vs.\nFallow", "vs. Fallow", c_song$Comparisonlab)
c_song$Comparisonlab <- gsub("vs.\nContinuous", "vs. Continuous", c_song$Comparisonlab)
len <- dim(c_song)[1]
c_song$Comparison <- as.factor(c_song$Comparison)
c_song$Comparisonlab <- as.factor(c_song$Comparisonlab)
c_song$Comparison <- factor(c_song$Comparison, levels(c_song$Comparison)[c(2,4,1,3)])
c_song$Comparisonlab <- factor(c_song$Comparisonlab, levels(c_song$Comparisonlab)[c(2,4,1,3)])
p_song3 <- ggplot(dat=c_song, aes(x=Comparison, y=Est)) +
  scale_y_continuous(limits = ylim[1,]) +
  theme_bw() +
  theme(plot.margin = loc_left, panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  #geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="LRR", x="", title="G") +
  scale_x_discrete(labels =c_song$Comparisonlab) +
  geom_point(data=c_song[c(2:4),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est)) # refer to song 2018 for significance
p_song3




# strom carbon
c_strom$Comparisonlab <- gsub(" ", "\n", c_strom$Comparison)
c_strom <- c_strom[order(c_strom$Comparison, c_strom$Year),]
len <- dim(c_strom)[1]
c_strom$X <- c(1,1,2,2)
p_strom2 <- ggplot(dat=c_strom, aes(x=Comparison, y=Est)) +
  scale_y_continuous(limits = ylim[1,], sec.axis = second.axis) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal", plot.margin = loc_right,  axis.text.y.left = element_blank(), panel.grid.major = element_line(colour=gridcol), panel.border=element_rect(colour=boxcol, fill=NA, size=boxthick), strip.text.x = element_text(size = stripsize), axis.title=element_text(size=axissize), axis.text=element_text(size=axistext)) +
  geom_hline(yintercept=0, linetype=zerolinetype, color=zerolinecol, size=zerolinethick) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), colour="black", width=0, size=1) +
  geom_point(size=ptsize, shape=ptshape, fill=ptbg) +
  labs(y="", x="", title="H") +
  scale_x_discrete(labels =c_strom$Comparisonlab[c(1,3)]) +
  geom_point(data=c_strom[which(sign(c_strom$CI_lower)==sign(c_strom$CI_upper)),], shape=ptshapesig, size=ptsize, aes(x=X, y=Est)) +
  facet_grid(.~Year, labeller = label_both) 
p_strom2





# Plot titles
tash_m <- ggplot() + lims(x=c(-10,10), y=c(-2,2)) +
  annotate(geom="text", y=1, x=1, label=expression(paste("Ashworth ",italic("et al.")," (2017)")), size=txsize, angle = 0) + 
  annotate(geom="text", y=-1, x=1, label=expression(paste("MTREC")), size=txsize, angle = 0) + theme_void()
tash_r <- ggplot() +  lims(x=c(-10,10), y=c(-2,2)) +
  annotate(geom="text", y=1, x=-1, label=expression(paste("Ashworth ",italic("et al.")," (2017)")), size=txsize, angle = 0) + 
  annotate(geom="text", y=-1, x=-1, label=expression(paste("RECM")), size=txsize, angle = 0) + theme_void()
tcloutier <- ggplot() + lims(x=c(-10,10), y=c(-2,2)) +
  annotate(geom="text", x=1, y=0, label=expression(paste("Cloutier ",italic("et al.")," (2020)")), size=txsize, angle = 0) + theme_void()
tgao <- ggplot() + lims(x=c(-10,10), y=c(-2,2)) +
  annotate(geom="text", x=-1, y=0, label=expression(paste("Gao ",italic("et al.")," (2019)")), size=txsize, angle = 0) + theme_void()
tsong <- ggplot() + lims(x=c(-10,10), y=c(-2,2)) +
  annotate(geom="text", x=1, y=0, label=expression(paste("Song ",italic("et al.")," (2018)")), size=txsize, angle = 0) + theme_void()
tstrom <- ggplot() + lims(x=c(-10,10), y=c(-2,2)) +
  annotate(geom="text", x=-1, y=0, label=expression(paste("Strom ",italic("et al.")," (2020)")), size=txsize, angle = 0) + theme_void()

# NA
plot_NA <- ggplot() + lims(y=c(-10,10), x=c(-2,2)) +
  annotate(geom="text", x=0, y=0, label="NA", size=5) + theme_void()

# plot nothing
plot_no <- ggplot() + lims(y=c(-10,10), x=c(-2,2)) + theme_void()






pdf("Figures/effect-plots/*1_main-fig_LRR-carbon.pdf", height=11, width=11)
laymat <- rbind(c(1,1,2,2),
                c(3,4,5,6),
                c(7,7,8,8),
                c(9,9,10,10),
                c(11,11,12,12),
                c(13,13,14,14))
gridExtra::grid.arrange(tash_m, tash_r, 
                        p_ashworth3, p_ashworth4, p_ashworth7, p_ashworth8, 
                        tcloutier, tgao,
                        p_cloutier2, p_gao3, 
                        tsong, tstrom,
                        p_song3, p_strom2, 
                  nrow=6, heights=c(0.15,1,0.15,1,0.15,1),
                  ncol=4, widths=c(1.1,0.9,0.9,1.1),
                  layout_matrix=laymat)
dev.off()
 





