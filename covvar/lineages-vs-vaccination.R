rm(list = ls())

# version 3.5.0
# Reference: Wuhan-Hu-1 with XBB SNPs
# Updated at: 2024-04-25 01:03:07 (UTC)
# Dataset name: nextstrain/sars-cov-2/XBB

setwd("~/Desktop/COVAR/pangolin")

nextclade <- read.csv("~/Desktop/COVAR/pangolin/nextclade-XBB2.csv", sep=";")

nextclade <- nextclade[c('seqName', 'Nextclade_pango',
                         'qc.overallStatus', 'totalSubstitutions',
                         'totalAminoacidSubstitutions',
                         'aaSubstitutions')]

dates<-read.csv('covvar-dates_updated.csv')
head(dates)

nextclade$seqName<-gsub('_consensus_NC_045512_2', '', nextclade$seqName)
nextclade$seqName<-gsub('Project_PE_', '', nextclade$seqName)

nextclade$seqName<-gsub('-', '_', nextclade$seqName)
nextclade$seqName <- gsub('_S20_L001_001 posreads: 10; cutoff: 0.1; mixrate: 0.0005', '', nextclade$seqName)
nextclade$seqName<-gsub('C24_01_007_', '', nextclade$seqName)

dates$COVVAR_ID<-gsub('-', '_', dates$COVVAR_ID)

nextclade$Month <- dates$Dates[match(nextclade$seqName, dates$COVVAR_ID)]
nextclade$Year <- dates$Year[match(nextclade$seqName, dates$COVVAR_ID)]

patterns <- c('XBB.1', 'XBB.2', 'JN.1', 'FL', 'EG', 'CM')

nextclade$lineage <- ifelse(
  grepl(patterns[1], nextclade$Nextclade_pango), patterns[1],
  ifelse(
    grepl(patterns[2], nextclade$Nextclade_pango), patterns[2],
    ifelse(
      grepl(patterns[3], nextclade$Nextclade_pango), patterns[3],
      ifelse(
        grepl(patterns[4], nextclade$Nextclade_pango), patterns[4],
        ifelse(
          grepl(patterns[5], nextclade$Nextclade_pango), patterns[5],
      ifelse(
        grepl(patterns[6], nextclade$Nextclade_pango), patterns[6],
        "Others"
      )
    )
  )
)
)
)
table(nextclade$lineage)
# CM     EG     FL   JN.1 Others  XBB.1  XBB.2 
# 4      7     11      4      5    113     11

nextclade$date <- paste0(nextclade$Month, " ", nextclade$Year)
nextclade$date <- as.Date(paste("01", nextclade$date), format = "%d %b %Y")
nextclade$date2 <- format(nextclade$date, "%b %y")
nextclade$date2

nextclade_short <- nextclade

lineage_df <- data.frame(table(nextclade$lineage, nextclade$date2))
lineage_df <- lineage_df[lineage_df$Freq>0,]
head(lineage_df)

tmp_df <- data.frame(Var1=c('FL','Others'),
                     Var2=c('Aug 23', 'Sep 23'),
                     Freq=c(0,0))
lineage_df <- rbind(lineage_df, tmp_df)
monthsyr <- c('Apr 23', 'May 23', 'Jun 23', 'Jul 23', 'Aug 23', 'Sep 23', 'Oct 23', 'Dec 23', 'Jan 24', 'Feb 24')

lineage_df$Var2 <- factor(lineage_df$Var2,
                          levels = monthsyr,
                          labels = monthsyr)

lineage_df$Var1 <- factor(lineage_df$Var1,
                          levels = c('FL','Others', 'CM', 'XBB.1', 'XBB.2', 'EG', 'JN.1')
                          )

p<-ggplot(lineage_df, aes(x=Var2, y=Freq,fill=Var1)) + #fill="#4682B4"
 # facet_wrap(~Var2, scales = 'free_x')+
  geom_bar(stat = "identity")+ylab("Number of sequences")+xlab("")+
 # geom_text(aes(label=paste0(Freq+4)))+
  theme(axis.title.x = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=12))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Status",values=rev(c("#CD7F32","#4682B4"))) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.2, hjust = 0.95),
        panel.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=12),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12),
        legend.position = "right"
  )+
  scale_y_continuous(labels = function(x) format(x, digits=1))+
scale_fill_brewer(palette = "Set3", direction = 1)#+
#ggtitle(label = "Virus family: Filoviridae")+theme(legend.position = "none")
p$labels$fill<-"Lineages"
p

muts_df <- nextclade[c('seqName', 'lineage', 'totalAminoacidSubstitutions', 'totalSubstitutions')]
names(muts_df)<-c('seqName', 'lineage', 'AA mutations', 'NT mutations')
names(muts_df)


muts_df$lineage <- factor(muts_df$lineage,
                          levels = c('FL','Others', 'CM', 'XBB.1', 'XBB.2', 'EG', 'JN.1'))

long_data <- reshape2::melt(muts_df, id.vars = c("seqName","lineage"),
                            measure.vars=c("AA mutations",
                                                 'NT mutations'), variable.name = "Variable", value.name = "Value")


p1<-ggplot(long_data, aes(x=lineage, y=Value,fill=lineage)) + #fill="#4682B4"
  facet_wrap(~Variable, scales = 'free_x', nrow = 2)+ylab("Number of mutations")+xlab("")+
  geom_boxplot(outlier.shape = NA)+
  # geom_text(aes(label=paste0(Freq+4)))+
  theme(axis.title.x = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=12))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Status",values=rev(c("#CD7F32","#4682B4"))) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.2, hjust = 0.95),
        panel.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=12),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12),
        legend.position = "right"
  )+
  scale_y_continuous(labels = function(x) format(x, digits=1))+
  scale_fill_brewer(palette = "Set3", direction = 1)#+
#ggtitle(label = "Virus family: Filoviridae")+theme(legend.position = "none")
p1$labels$fill<-"Lineages"
p1

aamuts <- nextclade[c('seqName', 'lineage','aaSubstitutions')]

library(tidyr)
aamuts_long <- aamuts %>%
  separate_rows(aaSubstitutions, sep = ",")

aamuts_long$count <- 1

df3 <- aggregate(count~aaSubstitutions, data = aamuts_long, FUN = sum)
df3$protein <- sub(":.*", "", df3$aaSubstitutions)
df3$mutation <- sub(".*:", "", df3$aaSubstitutions)
table(df3$protein)

df4 <- aggregate(count~aaSubstitutions+lineage, data = aamuts_long, FUN = sum)
df4$protein <- sub(":.*", "", df4$aaSubstitutions)
df4$mutation <- sub(".*:", "", df4$aaSubstitutions)
table(df4$protein)

spike <- df4[df4$protein=="S",]
spike2 <- aggregate(count~mutation+lineage, data = spike, FUN = sum)

names(spike)

unique_mutations <- spike %>%
  group_by(lineage) %>%
  filter(!mutation %in% unlist(spike[-which(spike$lineage == unique(lineage)), "mutation"]))


unique_mutations$lineage <- factor(unique_mutations$lineage,
                          levels = c('FL','Others', 'CM', 'XBB.1', 'XBB.2', 'EG', 'JN.1'))

# p2<-ggplot(unique_mutations, aes(x=aaSubstitutions, y=count,fill=lineage)) + #fill="#4682B4"
#   facet_wrap(~lineage, scales = 'free_y', nrow = 1)+
#   geom_bar(stat = "identity")+ylab("Number of sequences")+xlab("")+
#   theme(axis.title.x = element_text(color="black", size=15))+
#   theme(axis.text.x = element_text(color="black", size=12))+
#   theme(axis.title.y = element_text(color="black", size=15))+
#   theme(axis.text.y = element_text(color="black", size=12))+
#   #scale_fill_manual("Status",values=rev(c("#CD7F32","#4682B4"))) +
#   theme(axis.text.x = element_text(angle=90, vjust = 0.2, hjust = 0.95),
#         panel.background = element_blank(),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=12),
#         strip.background = element_rect(fill = "white", colour = "black"),
#         strip.text = element_text(size=15),#
#         panel.border = element_rect(fill=NA),
#         title = element_text(size = 12),
#         legend.position = "none"
#   )+
#   scale_y_continuous(labels = function(x) format(x, digits=1))+
#   scale_fill_brewer(palette = "Set3", direction = 1)#+
# #ggtitle(label = "Virus family: Filoviridae")+theme(legend.position = "none")
# p2$labels$fill<-"Lineages"
# p2

aamuts_long$protein <- sub(":.*", "", aamuts_long$aaSubstitutions)
aamuts_long$mutation <- sub(".*:", "", aamuts_long$aaSubstitutions)
aamuts_long_spike <- aamuts_long[aamuts_long$protein=="S",]
aamuts_long_spike <- aamuts_long_spike[c('seqName','lineage', 'mutation', 'count')]
aamuts_long_spike<-aamuts_long_spike[!duplicated(aamuts_long_spike),]
df5 <- aggregate(count~seqName+lineage, data = aamuts_long_spike, FUN = sum)

df5$lineage <- factor(df5$lineage,levels = c('FL','Others', 'CM', 'XBB.1', 'XBB.2', 'EG', 'JN.1'))
df5$label <- 'Spike AA mutations'
p3<-ggplot(df5, aes(x=lineage, y=count,fill=lineage)) + #fill="#4682B4"
  facet_wrap(~label)+
  geom_boxplot(outlier.shape = NA)+ylab("Number of mutations")+xlab("")+
  theme(axis.title.x = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=12))+
  theme(axis.text.y = element_text(color="black", size=12))+
  #scale_fill_manual("Status",values=rev(c("#CD7F32","#4682B4"))) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.2, hjust = 0.95),
        panel.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=12),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12),
        legend.position = "right"
  )+
  scale_y_continuous(labels = function(x) format(x, digits=1))+
  scale_fill_brewer(palette = "Set3", direction = 1)#+
#ggtitle(label = "Virus family: Filoviridae")+theme(legend.position = "none")
p3$labels$fill<-"Lineages"
p3


img1 <- cowplot::plot_grid(p,p3, nrow = 2)
img1

cowplot::plot_grid(img1, p1, nrow = 1)

### Visualising the phylogenetic tree

require(ggtreeExtra)
treef <- read.tree('/Users/alfred/Desktop/COVAR/pangolin/run002/clean-aln.nwk')
str(treef)

rooted_tree <- midpoint.root(treef)

p <- ggtree(rooted_tree)+#geom_tiplab(size=3)+
 ylim(-50, 320) + theme_tree2()
p
p<-p+geom_segment(aes(x = 0.001, y = -40, xend = 0.002, yend = -40))
p<-p+geom_text(aes(x = 0.0015, y = -30, label='0.001'))

  
nextclade <- read.csv2("~/Desktop/COVAR/fasttree/nxt2.tsv", sep = '\t')
head(nextclade)
#nextclade <- nextclade[c(1,3)]
rownames(nextclade) <- nextclade$seqName
nextclade$seqName<-NULL  

nextclade$group2<-nextclade$group
nextclade$group[nextclade$group2=='hCoV-19']<-"GISAID"
nextclade$group[nextclade$group2!='hCoV-19']<-"COVVAR"

library(ggnewscale)
names(nextclade)<-c("Lineages", "Source", "Group")

p2<-gheatmap(p, nextclade[1],
         colnames=TRUE,
         colnames_position = "bottom",
         offset = 0.0001,
         width=0.2,
         font.size = 4,
         colnames_angle=-45, hjust = 0,
         colnames_offset_y = -5) +
  scale_fill_brewer(palette = 'Accent', direction = 1) 
p2$labels$fill<-'Lineages'

p3 <- p2 + new_scale_fill()
p4<-gheatmap(p3, nextclade[2],
         colnames=TRUE,
         colnames_position = "bottom",
         offset = 0.0008,
         width=0.2,
         font.size = 4,
         colnames_angle=-45, hjust = 0,
         colnames_offset_y = -5) +
  scale_fill_brewer(palette = 'Pastel1', direction = 1) 
p4$labels$fill<-"Source"
p4 <- p4+theme_void()
p4

######### vaccination vs lineages

nextclade_short$seqName<-gsub('_S20','',nextclade_short$seqName)
nextclade_short$seqName<-gsub('_RPT','',nextclade_short$seqName)

vac_status <- read.csv("WHOAFROMoVECOVID_DATA_2024-04-12_0922_KR.csv")
vac_status <- vac_status[c('study_id','cov_vacination')]

vac_status$study_id <- gsub('-', '_', vac_status$study_id)

vac_status <- vac_status[vac_status$study_id%in%nextclade_short$seqName, ]

vac_status$lineage <- nextclade_short$lineage[match(
  vac_status$study_id, nextclade_short$seqName
)]
table(vac_status$cov_vacination)

vac <- data.frame(table(vac_status$cov_vacination))
vac
vac$percent <- vac$Freq / sum(vac$Freq) * 100
vac<-vac[order(vac$Freq, decreasing = TRUE),]
vac$Var1 <- factor(vac$Var1, levels = c(1,2,3), 
                   labels =c('Vaccinated', 'Not vaccinated', 'Unknown'))

vacp1 <- ggplot2::ggplot(data = vac, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +# ylim(0,120)+
  geom_text(aes(y=Freq+3, label = paste0(Freq, " (", round(percent), "%)"))) +
  ylab("Number of participants")+xlab("Vaccination status")+
  theme(axis.title.x = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=12))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(#axis.text.x = element_text(angle=0, vjust = 0.2, hjust = 0.95),
    panel.background = element_blank(),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size=12),#
    panel.border = element_rect(fill=NA),
    title = element_text(size = 12),
    legend.position = "none"
  )+#coord_flip()+
  scale_fill_brewer(palette = "Set2", direction = -1)
vacp1

vac_status$study_id[!vac_status$study_id%in%nextclade_short$seqName]
nextclade_short$seqName[!nextclade_short$seqName%in%vac_status$study_id]

#vac_status <- vac_status[vac_status$cov_vacination_status%in%c(1,2),]
vac_lineages<-table(vac_status$lineage, vac_status$cov_vacination)
vac_lineages

chisq.test(vac_lineages)
fisher.test(vac_lineages)

vac_lineages<-data.frame(vac_lineages)
names(vac_lineages)<-c('Lineage','Status','Freq')
vac_lineages <- vac_lineages[vac_lineages$Freq>0,]
vac_lineages$Status <- factor(vac_lineages$Status, levels = c(1,2,3), 
                   labels =c('Vaccinated', 'Not vaccinated', 'Unknown'))
vacp2 <- ggplot2::ggplot(data = vac_lineages, aes(x = Lineage, y = Freq, fill=Status)) +
  geom_bar(stat = "identity") +
  ylab("Number of participants")+xlab("Lineages")+
  theme(axis.title.x = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=12))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(#axis.text.x = element_text(angle=0, vjust = 0.2, hjust = 0.95),
    panel.background = element_blank(),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size=12),#
    panel.border = element_rect(fill=NA),
    title = element_text(size = 12),
    legend.position = "right"
  )+#coord_flip()+
  scale_fill_brewer(palette = "Set2", direction = -1)
vacp2

img2 <- cowplot::plot_grid(vacp1, vacp2, nrow = 1, rel_widths = c(1,1.5))
img2

wide_data <- reshape2::dcast(vac_lineages, Lineage ~ Status, value.var = "Freq")
wide_data
#write.csv(wide_data, 'vaccine-status-lineages.csv', quote = FALSE, row.names = FALSE)

#CM     EG     FL   JN.1 Others  XBB.1  XBB.2 

cm<-nextclade_short$Nextclade_pango[nextclade_short$lineage=='CM']
table(cm)
eg<-nextclade_short$Nextclade_pango[nextclade_short$lineage=='EG']
table(eg)
fl<-nextclade_short$Nextclade_pango[nextclade_short$lineage=='FL']
table(fl)
jn<-nextclade_short$Nextclade_pango[nextclade_short$lineage=='JN.1']
table(jn)
others<-nextclade_short$Nextclade_pango[nextclade_short$lineage=='Others']
table(others)
xbb1<-nextclade_short$Nextclade_pango[nextclade_short$lineage=='XBB.1']
table(xbb1)
xbb2<-nextclade_short$Nextclade_pango[nextclade_short$lineage=='XBB.2']
table(xbb2)

