#Read the configuration file. Please avoid numeric labels in the alignment!

library(ggplot2)
barcoding.gap.cfg <- scan(file="../barcoding_gap.cfg",what="character",strip.white=TRUE,quiet=TRUE)

#Determines the number of groups, the file to be read, and the number of taxa.

group.names <- character()
for (i in 1:length(barcoding.gap.cfg)) {
	ifelse(barcoding.gap.cfg[i] == "###Group###",group.names <- c(group.names,barcoding.gap.cfg[i+1]),NA)
	ifelse(barcoding.gap.cfg[i] == "###Filename###",file.name <- barcoding.gap.cfg[i+1],NA)
	}
num.groups <- length(group.names)
group.names.string <- group.names[1]
for (i in 2:num.groups) {
	group.names.string <- paste(group.names.string,group.names[i],sep="_")
	}

#Read the distance file outputted by distmat and builds a matrix using "NA" for blanks. Change 8 with 9 in the following line if you used the -ambiguous option of distmat.

barcoding.gap.distmat <- scan(file=paste("../",file.name,".distmat",sep=""),skip=8,what="character",strip.white=TRUE,quiet=TRUE)
num.taxa <- (-5+sqrt(25+8*length(barcoding.gap.distmat)))/2
barcoding.gap.distances <- numeric()
to.be.eliminated <- numeric()
current.last <- -2
for (i in 1:num.taxa) {
	current.first <- current.last+3
	current.last <- current.first+num.taxa-i
	barcoding.gap.distances <- c(barcoding.gap.distances,rep(NA,i-1))
	barcoding.gap.distances <- c(barcoding.gap.distances,as.numeric(barcoding.gap.distmat[current.first:current.last]))
	ifelse(is.element(barcoding.gap.distmat[current.last+1],barcoding.gap.cfg),NA,to.be.eliminated <- c(to.be.eliminated,i))
	}
barcoding.gap.matrix <- matrix(barcoding.gap.distances,ncol=num.taxa,nrow=num.taxa,byrow=TRUE)
barcoding.gap.matrix.backup <- barcoding.gap.matrix
for (i in 1:length(to.be.eliminated)) {
	barcoding.gap.matrix[to.be.eliminated[i],] <- NA
	barcoding.gap.matrix[,to.be.eliminated[i]] <- NA
	}

#Determines numbers of taxa of each group and discriminates intra-group and inter-group distances.

correspondences.taxon <- character()
correspondences.number <- numeric()
intra.group <- numeric()
inter.group <- numeric()
first.taxon <- numeric()
second.taxon <- numeric()
starts.ends <- numeric()
for (i in 1:length(barcoding.gap.cfg)) {
	ifelse(barcoding.gap.cfg[i] == "###Group###",starts.ends <- c(starts.ends,i-1,i+2),NA)
	}
starts.ends <- c(starts.ends[2:length(starts.ends)],length(barcoding.gap.cfg))
for (i in 1:num.groups) {
	current.group.taxa <- barcoding.gap.cfg[starts.ends[i*2-1]:starts.ends[i*2]]
	current.group.pos <- numeric()
	for (j in 1:length(barcoding.gap.distmat)) {
		if (is.element(barcoding.gap.distmat[j],current.group.taxa)) {
			current.group.pos <- c(current.group.pos,as.numeric(barcoding.gap.distmat[j+1]))
			correspondences.taxon <- c(correspondences.taxon,barcoding.gap.distmat[j])
			correspondences.number <- c(correspondences.number,as.numeric(barcoding.gap.distmat[j+1]))
			}
		}
	for (j in 1:length(current.group.pos)) {
		for (k in j:length(current.group.pos)) {
			if (j == k) NA
				else {
				intra.group <- c(intra.group,barcoding.gap.matrix[current.group.pos[j],current.group.pos[k]])
				first.taxon <- c(first.taxon,current.group.pos[j])
				second.taxon <- c(second.taxon,current.group.pos[k])
				}
			barcoding.gap.matrix[current.group.pos[j],current.group.pos[k]] <- NA
			}
		}
	}
for (j in 1:num.taxa) {
	for (k in j:num.taxa) {
		if (is.na(barcoding.gap.matrix[j,k]) || j == k) NA
			else {
			inter.group <- c(inter.group,barcoding.gap.matrix[j,k])
			first.taxon <- c(first.taxon,j)
			second.taxon <- c(second.taxon,k)
			}
		}
	}

#Finally, the histogram. Be brave.

correspondences.data.frame <- data.frame(number=correspondences.number,taxon=correspondences.taxon)
correspondences.data.frame <- correspondences.data.frame[order(correspondences.data.frame$number),]
distance.data.frame <- data.frame(first=first.taxon,second=second.taxon,dis.type=c(rep("intra",length(intra.group)),rep("inter",length(inter.group))),distance=c(intra.group,inter.group))
barcoding.gap.histogram <- ggplot(distance.data.frame,aes(distance,fill=dis.type))
barcoding.gap.histogram <- barcoding.gap.histogram + scale_fill_manual(values=c("orangered2","chartreuse4"))
barcoding.gap.histogram <- barcoding.gap.histogram + geom_density(alpha=0.8)
barcoding.gap.histogram <- barcoding.gap.histogram + xlim(0,ceiling(max(distance.data.frame$distance)))
barcoding.gap.histogram <- barcoding.gap.histogram + theme(line=element_line(colour="black",size=0.1),axis.line=element_line(colour="black",size=0.1),axis.ticks=element_line(colour="black",size=0.1),panel.grid=element_blank(),panel.background=element_rect(fill="white",colour="black"),plot.background=element_rect(fill="white",colour="black"),legend.background=element_blank(),legend.title=element_blank(),text=element_text(family="Times",colour="black"),axis.text=element_text(family="Times",colour="black"))
barcoding.gap.histogram <- barcoding.gap.histogram + ggtitle(paste("Dataset: ",file.name,"\nGroups: ",group.names.string,sep=""))
print(barcoding.gap.histogram)
dev.copy2pdf(file=paste("./barcoding_gap_",file.name,"_",group.names.string,".pdf",sep=""))
write.table(distance.data.frame,file=paste("./barcoding_gap_",file.name,"_",group.names.string,".txt",sep=""),append=FALSE,quote=FALSE,sep="\t",eol="\n",row.names=FALSE,col.names=FALSE)
write.table(correspondences.data.frame,file=paste("./correspondences_",file.name,"_",group.names.string,".txt",sep=""),append=FALSE,quote=FALSE,sep="\t",eol="\n",row.names=FALSE,col.names=TRUE)

#I still have to implement the computation of shadowed areas!...
