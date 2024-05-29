rm(list = ls())

require(ggplot2)
library(dplyr)
require(cowplot)

setwd("/path/to/COVAR")

# import nextclade results into the R environment
nextclade <- read.csv("nextclade-XBB.csv", sep=";")

# version 3.5.0
# Reference: Wuhan-Hu-1 with XBB SNPs
# Updated at: 2024-04-25 01:03:07 (UTC)
# Dataset name: nextstrain/sars-cov-2/XBB

# pick the relevant columns - we are interested in lineage assignment, QC and mutations both nt and aa 
nextclade <- nextclade[c('seqName', 'Nextclade_pango',
                         'qc.overallStatus', 'totalSubstitutions',
                         'totalAminoacidSubstitutions',
                         'aaSubstitutions')]

# import sample collection dates for the analysed sequences
dates<-read.csv('covvar-dates_updated.csv')
head(dates)

# sequence names are all in GISAID form, we need to pick jusy the sample ID so we are able to match it with the IDs in collection dates
nextclade <- nextclade %>%
  mutate(seqName = sub("hCoV-19/Uganda/", "", seqName)) %>%
  mutate(seqName = sub("/2023", "", seqName)) %>%
  mutate(seqName = sub("/2024", "", seqName)) %>%
  mutate(seqName = sub("/hCoV-19", "", seqName)) %>%
mutate(seqName = sub("/", "", seqName))

# we replace hyphen (-) with underscore (_) in seqnames and samples IDs in nexclade and date dataset respectively
nextclade$seqName<-gsub('-', '_', nextclade$seqName)
dates$COVVAR_ID<-gsub('-', '_', dates$COVVAR_ID)

# add sampling dates to nextclade dataset by matching IDs
nextclade$Month <- dates$Dates[match(nextclade$seqName, dates$COVVAR_ID)]
nextclade$Year <- dates$Year[match(nextclade$seqName, dates$COVVAR_ID)]

# create a vector of patterns to look up in the lineage assignment
patterns <- c('XBB.1', 'XBB.2', 'JN.1', 'FL', 'EG', 'CM')

# create a new column with new lineages depending on patterns above - e.g XBB.1 will be for all XBB.1 and respective sublineages
# everything else that is not part of the above patterns will be classified as Others 

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

# do some date formating to enable us visualise the dates well
nextclade$date <- paste0(nextclade$Month, " ", nextclade$Year)
nextclade$date <- as.Date(paste("01", nextclade$date), format = "%d %b %Y")
nextclade$date2 <- format(nextclade$date, "%b %y")
nextclade$date2

# create a copy of this nexclade dataset before further manipulation - we shall refer to this later while working on the vaccination vs lineages
nextclade_short <- nextclade

# create a dataframe of number of sequences in a particular lineage across the sampling dates
lineage_df <- data.frame(table(nextclade$lineage, nextclade$date2))
lineage_df <- lineage_df[lineage_df$Freq>0,]
head(lineage_df)

# we did realise we had no data in Aug 23 and Sep 23, as such we add some entries manually with zero counts so that we can show this on the plots 
tmp_df <- data.frame(Var1=c('FL','Others'),
                     Var2=c('Aug 23', 'Sep 23'),
                     Freq=c(0,0))
lineage_df <- rbind(lineage_df, tmp_df)

# we create a vector in the order of dates to be used for visualisation
monthsyr <- c('Apr 23', 'May 23', 'Jun 23', 'Jul 23', 'Aug 23', 'Sep 23', 'Oct 23', 'Dec 23', 'Jan 24', 'Feb 24', 'Mar 24')
lineage_df$Var2 <- factor(lineage_df$Var2,
                          levels = monthsyr,
                          labels = monthsyr)

# similarly we create the order in which to dispaly the lineages
lineage_df$Var1 <- factor(lineage_df$Var1, levels = c('FL','Others', 'CM', 'XBB.1', 'XBB.2', 'EG', 'JN.1'))

# plot the lineages vs sampling dates
p<-ggplot(lineage_df, aes(x=Var2, y=Freq,fill=Var1)) + 
  geom_bar(stat = "identity")+ylab("Number of sequences")+xlab("")+
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
scale_fill_brewer(palette = "Set3", direction = 1)
p$labels$fill<-"Lineages"
p

# now we pick the lineages and corresponsing number of nt and aa substitution mutations for all sequences analysed
muts_df <- nextclade[c('seqName', 'lineage', 'totalAminoacidSubstitutions', 'totalSubstitutions')]
names(muts_df)<-c('seqName', 'lineage', 'AA mutations', 'NT mutations')

# similarly we create the order in which to dispaly the lineages
muts_df$lineage <- factor(muts_df$lineage, levels = c('FL','Others', 'CM', 'XBB.1', 'XBB.2', 'EG', 'JN.1'))

# prepare data for plotting in ggplot2 - basically creating a long format representation of the above data
long_data <- reshape2::melt(muts_df, id.vars = c("seqName","lineage"),measure.vars=c("AA mutations",'NT mutations'), variable.name = "Variable", value.name = "Value")

# plot the distribution of the number of aa and nt mutations across lineages                  
p1<-ggplot(long_data, aes(x=lineage, y=Value,fill=lineage)) + 
  facet_wrap(~Variable, scales = 'free_x', nrow = 1)+ylab("Number of mutations")+xlab("")+
  geom_boxplot(outlier.shape = NA)+
  theme(axis.title.x = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=12))+
  theme(axis.text.y = element_text(color="black", size=12))+
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
p1$labels$fill<-"Lineages"
p1

# over to the aa mutations                     
aamuts <- nextclade[c('seqName', 'lineage','aaSubstitutions')]
# create a long form representation of these data                     
aamuts_long <- aamuts %>%
  separate_rows(aaSubstitutions, sep = ",")
# add a column to specify the protein and mutations respectively
aamuts_long$protein <- sub(":.*", "", aamuts_long$aaSubstitutions)
aamuts_long$mutation <- sub(".*:", "", aamuts_long$aaSubstitutions)
# pick only data for the spike protein                     
aamuts_long_spike <- aamuts_long[aamuts_long$protein=="S",]
# add a pseudo counter on this dataframe 
aamuts_long_spike$count <- 1               
aamuts_long_spike <- aamuts_long_spike[c('seqName','lineage', 'mutation', 'count')]
# remove any duplicates before aggregation                     
aamuts_long_spike<-aamuts_long_spike[!duplicated(aamuts_long_spike),]
# compute the number of spike mutations in each sequence - do this per lineage
df5 <- aggregate(count~seqName+lineage, data = aamuts_long_spike, FUN = sum)
# specify the order in which lineages should be displayed on the plot
df5$lineage <- factor(df5$lineage,levels = c('FL','Others', 'CM', 'XBB.1', 'XBB.2', 'EG', 'JN.1'))
df5$label <- 'Spike AA mutations'
# plot the distribution of spike aa mutations across lineages                    
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

# group images that we have generated so far into a single plot
img1 <- cowplot::plot_grid(p,p3, nrow = 2)
cowplot::plot_grid(img1, p1, nrow = 1)

######### vaccination vs lineages ##########
                     
# import data on vaccination status  
vac_status <- read.csv("WHOAFROMoVECOVID_DATA_2024-04-12_0922_KR.csv")
# pick the relevant column - showing whether vaccinated or not                     
vac_status <- vac_status[c('study_id','cov_vacination')]
# repalce hypens with underscores in the study ID so theat they match with the nexclade IDs
vac_status$study_id <- gsub('-', '_', vac_status$study_id)
# pick vaccination data for which we have corresponding sequence data
vac_status <- vac_status[vac_status$study_id%in%nextclade_short$seqName, ]
# add lineage data to vacccine data
vac_status$lineage <- nextclade_short$lineage[match(
  vac_status$study_id, nextclade_short$seqName
)]
table(vac_status$cov_vacination)

# create a dataframe for vaccination status                     
vac <- data.frame(table(vac_status$cov_vacination))
vac$percent <- vac$Freq / sum(vac$Freq) * 100
vac<-vac[order(vac$Freq, decreasing = TRUE),]
vac$Var1 <- factor(vac$Var1, levels = c(1,2,3),labels =c('Vaccinated', 'Not vaccinated', 'Unknown'))
# plot the vaccination status data
vacp1 <- ggplot2::ggplot(data = vac, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +# ylim(0,120)+
  geom_text(aes(y=Freq+3, label = paste0(Freq, " (", round(percent,1), "%)"))) +
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

# create a contigency table showing vaccination status by lineage assignment                     
vac_lineages<-table(vac_status$lineage, vac_status$cov_vacination)
vac_lineages
# use Fisher's exact test to assess the association of vaccination status and lineage assignment
fisher.test(vac_lineages)
                     
# create dataframe to visualise relation beteen lineages and vaccination status
vac_lineages<-data.frame(vac_lineages)
names(vac_lineages)<-c('Lineage','Status','Freq')
vac_lineages <- vac_lineages[vac_lineages$Freq>0,]
vac_lineages$Status <- factor(vac_lineages$Status, levels = c(1,2,3),labels =c('Vaccinated', 'Not vaccinated', 'Unknown'))
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
# combine the two plots created around vaccine status and lineages 
img2 <- cowplot::plot_grid(vacp1, vacp2, nrow = 1, rel_widths = c(1,1.5))
img2
# create the data with counts of sequences - lineages by vaccination status
wide_data <- reshape2::dcast(vac_lineages, Lineage ~ Status, value.var = "Freq")
wide_data
#write.csv(wide_data, 'vaccine-status-lineages.csv', quote = FALSE, row.names = FALSE)

### Visualising the phylogenetic tree

# load the required packages                     
require(ggtreeExtra)
require(ggtree)
require(ape)
require(phytools)
require(ggnewscale)
                     
# import the treefile generated with IQtree
treef <- ape::read.tree('all-aln.fa.treefile')
str(treef)
# midpoint root the tree
rooted_tree <- midpoint.root(treef)
# create the first view of the tree
p <- ggtree(rooted_tree)+#geom_tiplab(size=3)+
 ylim(-50, 320) + theme_tree2()
p
# add the tree scale                      
p<-p+geom_segment(aes(x = 0.001, y = -40, xend = 0.002, yend = -40))
p<-p+geom_text(aes(x = 0.0015, y = -30, label='0.001'))

# load the nextclade results of the dataset used to generate the tree - this is a combination of background and study generated sequences
nextclade <- read.csv("nextclade-gisaid-cvr.csv", sep=";")
head(nextclade)
# of interest is just the lineage and the sequence name, so we pick that
nextclade<-nextclade[c(2,4)]
# we need to highlight study sequences (COVVAR) from background (GISAID)                     
nextclade$group[grepl('CVR',nextclade$seqName)]<-"COVVAR"
nextclade$group[is.na(nextclade$group)]<-"GISAID"
# create lineage patterns to higlight 
patterns <- c('BA', 'BQ.1', 'JN.1', 'XBB')
# assign a new lineage having specified patterns and associated sublineages
nextclade$lineage <- ifelse(
  grepl(patterns[1], nextclade$Nextclade_pango), patterns[1],
  ifelse(
    grepl(patterns[2], nextclade$Nextclade_pango), patterns[2],
    ifelse(
      grepl(patterns[3], nextclade$Nextclade_pango), patterns[3],
      ifelse(
        grepl(patterns[4], nextclade$Nextclade_pango), patterns[4],
        "Others"
      )
    )
  )
)
table(nextclade$lineage) # look at the lineage distribution

# we need row names of the annotation dataframe to be the same as the tip labels                     
rownames(nextclade) <- nextclade$seqName
nextclade$seqName<-NULL  
# now that we have a new assignment we remove the  Nextclade_pango column                    
nextclade$Nextclade_pango<-NULL
names(nextclade)<-c("Source", "Lineages")
# first we add the lineages annotation on the tree
p2<-gheatmap(p, nextclade[2],
         colnames=TRUE,
         colnames_position = "bottom",
         offset = 0.0001,
         width=0.05,
         font.size = 4,
         colnames_angle=-90, hjust = 0,
         colnames_offset_y = -5) +
  scale_fill_brewer(palette = 'Accent', direction = 1) 
p2$labels$fill<-'Lineages'
# second, we add the source (COVVAR or GISAID) annotation on the tree
p3 <- p2 + new_scale_fill()
p4<-gheatmap(p3, nextclade[1],
         colnames=TRUE,
         colnames_position = "bottom",
         offset = 0.00025,
         width=0.05,
         font.size = 4,
         colnames_angle=-90, hjust = 0,
         colnames_offset_y = -5) +
  scale_fill_brewer(palette = 'Pastel1', direction = 1) 
p4$labels$fill<-"Source"
p4 <- p4+theme_void()
p4

