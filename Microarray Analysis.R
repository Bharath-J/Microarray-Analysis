#Installing Bioconductor
#Bioconductor is a package that performs computational analysis 
# of biological problems

#This command runs an R program that defines the function
#biocLite that is used for installing Bioconductor
source("http://www.bioconductor.org/biocLite.R")

#Installing BioConductor packages "GEOquery" and "limma"
biocLite("GEOquery")
biocLite("limma")

#Activating the packages in the normal way
#Biobase is the core package of Bioconductor and is always 
#installed automatically
library(Biobase)
library(GEOquery)
library(limma)

#We also need gplots. Installing it in the normal way.
install.packages("gplots")

#Activating gplots
library(gplots)

#Now, lets download a dataset from NCBI Gene Expression Omnibus,
#http://www.ncbi.nlm.nih.gov/geo/
#The data for this task has the identifier "GDS1542" in this database.

#Downloading, putting it in the current directory and loading it into
#R in one step with the following command
gds <- getGEO('GDS1542', destdir=".")

#Let's see what's in the downloaded file.
show(gds)

#The command above shows data table with only 5 rows of data out of approximately
#22000 rows.

#Each row in the table are the measured value of one probe (spot)
#on the microarray.

#The cryptic name of this probe is in the first column (ID_REF).

#The gene (and protein produced from that gene) each probe detects 
#is in the second column.

#GEOData is a specialized class, for GEO derived datasets only.
#So, lets convert it into something that fits general functions of 
#Bioconductor.
eset <- GDS2eSet(gds)

#To get an over of eset
show(eset)

#Extracting the actual values as a standard matrix.
expdata <- exprs(eset)

#Filtering away the row/rows that has/have data missing.
w <- which(apply(is.na(expdata),1,sum)>0)
temp <- expdata[-w, ]
expdata <- temp

#Let's see how the data looks using a box plot
#But we should take 2-logarithm of data in order to make
#it look normally distributed
logdata <- log2(expdata)
boxplot(as.data.frame(logdata))
x11()

#To check whether the data is noisy when the signal is low,
#we will plot standard deviation against the mean for every probe.
probemeans <- apply(logdata, 1, mean)
probesd <- apply(logdata, 1, sd)
plot(probemeans, probesd)
x11()

#Eliminating probes that have weakest signals.
q25 <- quantile(probemeans, 0.25)
whichtosave <- which(probemeans > q25)
q25logdata <- logdata[whichtosave, ]

#Removing probes that have low variability.
mydata <- q25logdata[apply(q25logdata, 1, IQR)>1.5, ]

#Performing PCA - Principal Component Analysis on columns
#Taking transpose of the matrix
tdata <- aperm(mydata)
pca <- prcomp(tdata, scale=T)

#Plotting the samples in relation to the first two components
#and finding out which experimental condition they belong to
conditions <- phenoData(eset)$agent
plot(pca$x, type="n")
x11()
text(pca$x, labels=conditions, cex=0.5)

#Dendogram of correlation between the samples
pearsonCorr <- as.dist(1-cor(mydata))
hC <- hclust(pearsonCorr)
plot(hC, labels = conditions)
x11()

#Heatmap. Red -> high values, green -> low values
heatmap(mydata, col=greenred(100))


#Let's find the probes changing significantly
condfactor <- factor(eset$agent)

#Constructing a model matrix
design <- model.matrix(~0+condfactor)
colnames(design) <- c("ctrl", "tnf")

#Estimating the variances
fit <- lmFit(eset, design)

#Defining conditions we want to compare
contrastmatrix <- makeContrasts(tnf - ctrl, levels=design)

#Calculating p-values for the difference between the conditions
#defined by contrast matrix.
fit <- contrasts.fit(fit, contrastmatrix)
ebayes <- eBayes(fit)

#Histogram of the p-values
hist(ebayes$p.value)
x11()

#Regulating p-value cutoff in order to have control over
#false positives.
results <- decideTests(ebayes)

#Displaying the number of probes that passed
length(which(results!=0))

#Venn Diagram of results
vennDiagram(results)

#Extracting the original data for the most changing probes
resData <- exprs(eset)[results !=0, ]

#Adding gene symbols as row names
geneSymbol <- as.array(fData(eset)[,"Gene symbol"])
gs <- geneSymbol[c(which(results!=0))]
rownames(resData) <- gs

#Adding p-values in an extra column
pvalues <- ebayes$p.value[results!=0,]
resData <- cbind(resData, pvalues)

#Adding q-values
adj.pvalues <- p.adjust(ebayes$p.value, metho="BH")
adj.pvalues <- adj.pvalues[results!=0]
resData <- cbind(resData, adj.pvalues)

#Writing output to a file
write.table(resData, "most_regulated.txt", sep="\t", quote=FALSE)








