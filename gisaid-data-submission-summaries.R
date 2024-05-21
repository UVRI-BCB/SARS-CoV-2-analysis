#Script to generate summaries of data so far deposited in GISAID from Uganda
#Authors: Alfred Ssekagiri, Nicholas Bbosa
#Date created: 2021-09-28
#Date last modified: 2021-09-28

#Clear the environment
rm(list = ls())

#load required packages
require(ggplot2)
require(tidyverse)

#load the metadata obtained from GISAID, change the path to  a location on your machine
sarscov2<-read.delim("~/Desktop/gisaid_auspice_input_hcov-19_2021_09_28_10/1632825352353.metadata.tsv")
#retain the relevant columns for the summary (here we are interested in the lab, dates and pango-lineages)
sarscov2<-sarscov2[c("strain","date","originating_lab","pangolin_lineage")]

##################Figure 1: Timeline of sequencing efforts since beginning of outbreak ##############################
#create a new dataframe to be used for creating this figure, a copy the above sarcov2 dataset
df3<-sarscov2
#create a variable called number, set it to one, its is just a mere count of a sequence
df3$number<-1
#reinforce the date format
df3$Month_Yr_dd <- format(as.Date(df3$date), "%Y-%m-%d")
#create another date column with breaks of 1 week in between, call it date2
df3$date2<-as.Date(cut(as.Date(df3$Month_Yr_dd),
                            breaks = "week",
                            start.on.monday = FALSE))
#Now, obtain the genomes for which there sampling dates belong to period specified in the date2
df4<-aggregate(number~date2, data=df3, FUN = sum)
#Get the weekly cummulative sum of the genomes, store this in column cumm
df4$cumm<-cumsum(df4$number)
#Go on and generate the plot using ggplot2
sequencing_plot<-df4 %>%
  mutate(date=as.POSIXct(date2)) %>% #convert date to date
  ggplot(data=df4, mapping = aes(x = date2, y=cumm))+
  geom_line(width=20,color='blue')+geom_point(color='blue')+ #line plot with points on line, in blue
  theme(axis.text.x = element_text(color="black", size=16,angle = 90))+ #vertical alignment of dates and changing size
  scale_x_date(date_labels = "%d %b %Y", date_breaks = "4 week")+ #The date labels will have 4 weeks apart, this reduces cluter on the axis
  theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=15))+
  theme(legend.text = element_text(size=11))+
  theme(legend.title = element_text(size=12))+
  theme(legend.title = element_text(size=15, face="bold"))+
  theme(panel.background = element_blank(), #set a clear background with black border enclosing the figure
        panel.border = element_rect(fill=NA))+
  xlab('Sampling Date')+
  ylab('Number of Genomes')
sequencing_plot

############### Figure 2: Number of genomes by different laboratories #####################
#Here we are entirely working with the originating lab variable/column
#We make necessary changes as follows.
#Genomes from Rakai are essentially UVRI
sarscov2[sarscov2$originating_lab=="Rakai Health Sciences Program",]$originating_lab<-"UVRI"
#Genomes from Uganda Virus Research are labeled UVRI
sarscov2[sarscov2$originating_lab=="Uganda Virus Research Institute",]$originating_lab<-"UVRI"
#Genomes from Makerere are labeled MAK-CHS
sarscov2[sarscov2$originating_lab=="Integrated Biorepository of H3Africa Uganda – IBRH3AU",]$originating_lab<-"MAK-CHS"
#Genomes from CPHL are labeled UVRI
sarscov2[sarscov2$originating_lab=="Uganda Central Public Health Lab and Uganda Virus Research Institute",]$originating_lab<-"UVRI"
#Genomes from elsewhere other than those above is labeled MRC/UVRI and LSHTM
sarscov2[!sarscov2$originating_lab%in%c("CPHL & UVRI", "RHSP", "UVRI", "MAK-CHS"),]$originating_lab<-"MRC/UVRI & LSHTM"

#At this point we obtain the number of genomes by simply tabulating the 'originating_lab' column and making that a dataframe
df1<-data.frame(table(sarscov2$originating_lab))
#We add a column to the new dataframe to hold the percentage of genomes from each lab
df1$perc<-100*df1$Freq/sum(df1$Freq)

#get piechart of this information
pie = ggplot(df1, aes(x="", y=Freq, fill=Var1)) + geom_bar(stat="identity", width=1)#+facet_wrap(~Var1)
# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(Freq,"(",round(perc),"%",")")), 
                                                  position = position_stack(vjust = 0.5),
                                                  size=4)#, fontface="bold")
# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = "Var1")
# Tidy up the theme
pie = pie + theme_classic() + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
pie<-pie+theme(legend.text = element_text(size=12), #change font size of legend items here.
               legend.title = element_text(size = 15),#change font size of legend title here.
               title = element_text(hjust = 1.2)) #center align the title
pie$labels$fill<-"Laboratory"

pie

#get the equivalent barplot
barplot_genomes = ggplot(df1, aes(x=Var1, y=Freq, fill=Var1)) + geom_bar(stat="identity")#+facet_wrap(~Var1)
# Tidy up the theme
barplot_genomes<-barplot_genomes+theme(legend.text = element_text(size=12), #change font size of legend items here.
               legend.title = element_text(size = 15),
               axis.text = element_text(size=12),
               axis.title = element_text(size=15)) #center align the title
barplot_genomes<-barplot_genomes+ylab("Number of Genomes")+xlab("Laboratory")
barplot_genomes<-barplot_genomes+geom_text(aes(label = paste0(Freq,"(",round(perc),"%",")"),y=Freq+15), 
                              position = position_stack(),
                            size=4)
barplot_genomes<-barplot_genomes+theme(panel.background = element_blank(),
               panel.border = element_rect(fill=NA))

barplot_genomes$labels$fill<-"Laboratory"
barplot_genomes

############################Figure 3: Genomes and corresponding lineages/clades ##################
#Here we can inspect the pangolin-lineages column of sarscov2 dataframe
unique(sarscov2$pangolin_lineage)
#we would continue as get a data frame by tabulating as done before and categorise appropriately - Option 1.
#for this, Nicholas shared the breakdown following inspection and we did this manually (we were time bad) - Option 2.
#create a vector of lineage names
Var1<-c("Delta (B.1.6.17 & AY.4)", "A.23 & A.23.1","Eta", "Others")
#create a vector of corresponding number of genomes
Freq<-c(241, 233, 39,199)
#create a dataframe with two columns to accomodate the above information
df1<-data.frame(Var1=Var1, Freq=Freq)
#get the percentage
df1$perc<-100*df1$Freq/sum(df1$Freq)
#to maintain the same order, have the lineages as an ordered factor
df1$Var1 <- factor(df1$Var1, levels = df1$Var1)
#Generate the barplot
p = ggplot(df1, aes(x=Var1, y=Freq, fill=Var1)) + geom_bar(stat="identity")#+facet_wrap(~Var1)
# Tidy up the theme
p<-p+theme(legend.text = element_text(size=12), #change font size of legend items here.
               legend.title = element_text(size = 15),
               axis.text = element_text(size=12),
               axis.title = element_text(size=15)) #center align the title
p<-p+ylab("Number of Genomes")+xlab("")
p<-p+geom_text(aes(label = paste0(Freq,"(",round(perc),"%",")"),y=Freq+15), 
                   position = position_stack(),
                   size=4)
p<-p+theme(panel.background = element_blank(),panel.border = element_rect(fill=NA))
p$labels$fill<-"Lineage"
p
#get corresponding piechart
pie = ggplot(df1, aes(x="", y=Freq, fill=Var1)) + geom_bar(stat="identity", width=1)#+facet_wrap(~Var1)
# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(Freq,"(",round(perc),"%",")")), 
           position = position_stack(vjust = 0.5),
          size=4)#, fontface="bold")
# Add color scale (hex colors)
#pie = pie + scale_fill_manual(values=c("#55DDE0", "#33658A", "#2F4858", "#F6AE2D", "#F26419", "#999999")) 
# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = "Var1")
# Tidy up the theme
pie = pie + theme_classic() + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
pie<-pie+theme(legend.text = element_text(size=12), #change font size of legend items here.
               legend.title = element_text(size = 15),#change font size of legend title here.
               title = element_text(hjust = 1.2)) #center align the title
pie$labels$fill<-"Lineage"

pie


