args = commandArgs(trailingOnly=TRUE)
# Print the input file
message(paste("Input csv is:", args[1]))
message(paste("Output directory is:", args[3]))
message(paste("Append string is:", args[2]))
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).csv", call.=FALSE)
}
# Load data 
data <- read.csv(args[1], sep=",", skip = 1) 
# Libraries. Make sure they installed or install otherwise
if(!require(drc)){install.packages("drc", repos='http://cloud.r-project.org/' )}
if(!require(gridExtra)){install.packages("gridExtra", repos='http://cloud.r-project.org/')}
if(!require(grid)){install.packages("grid", repos='http://cloud.r-project.org/')}
if(!require(ggplot2)){install.packages("ggplot2", repos='http://cloud.r-project.org/')}
if(!require(gplots)){install.packages("gplots", repos='http://cloud.r-project.org/')}
if(!require(ggthemes)){install.packages("ggthemes", repos='http://cloud.r-project.org/')}
# Convert Well ID to character 
data$Well.ID <- as.character(data$Well.ID)
# I assume the format is cols 1-11 and cols 14-23 with 11 and 23 being the no antibody controls 
# I also assume that rows are paired BC, DE, etc .. 
# First row is B and last row is O
# Crete a matrix of row pairs
rows <- LETTERS[2:15]
pairs <- matrix(rows, ncol=2, byrow = T)
# Set the concentration vector in ng/ml
conc <- c(10000,2000, 400, 80, 16, 3.2, 0.64, 0.128, 0.0256)
conc.log <- log10(conc)
# Create list to store regression models and ggplots 
list.reg <- list()
list.rsd <- list()
list.eds <- list()
list.gfp.max <- list()
# Provide a vector of antibody names 
# For now its just AB1, AB2, etc ..
abs.left <- paste("AB_L",1:nrow(pairs), sep="_")
abs.right <- paste("AB_R",1:nrow(pairs), sep="_")
# Create color wheel 
col.wheel <- tableau_color_pal(palette = "tableau10")(10)
# 
message("Making left plate pdf")
# Start plotting pdfs
pdf(paste(args[3], args[2], "_left_plate_reg.pdf", sep=""),10,10)
# Iterate for every row pair 
for(i in 1:nrow(pairs)){
  # Subset data from main data frame 
  pair.1 <- data[grep(pairs[i,1], data$Well.ID),][1:10,] # 1:10 selects X02-X11
  pair.2 <- data[grep(pairs[i,2], data$Well.ID),][1:10,]
  # Get GFP value and normalize to no antibody control (last entry X11)
  # grep for 11 to get the X11 column from both replicates
  con.pair.1 <- pair.1$X.GFP.[grep(11, pair.1$Well.ID)] 
  con.pair.2 <- pair.2$X.GFP.[grep(11, pair.2$Well.ID)]
  # Correct the GFP values by dividing with the no antibody control 
  cor.pair.1.t <- pair.1$X.GFP./con.pair.1
  cor.pair.2.t <- pair.2$X.GFP./con.pair.2
  # Remove the zero value 
  cor.pair.1 <- cor.pair.1.t[-length(cor.pair.1.t)]
  cor.pair.2 <- cor.pair.2.t[-length(cor.pair.2.t)]
  # Get average values 
  cor.means <- rowMeans(cbind(cor.pair.1, cor.pair.2))
  cor.st <- apply(cbind(cor.pair.1, cor.pair.2), MARGIN = 1, range)
  # DOES IT MAKE SENSE ? 
  if(length(which(is.infinite(cor.means)==T)) ==0){
    # Fit using a 4-parameter log logistic model and store the model
    list.reg[[i]] <- drm(cor.means ~ conc,fct = LL.4())
    # Get ED5,10 and 50 values 
    list.eds[[i]] <- ED(object = list.reg[[i]], c(5,10,50), interval = "delta", display=F)
    # Also store %GFP and max antibody 
    list.gfp.max[[i]] <- cor.means[1]*100
    # Also store the residual standard error 
    list.rsd[[i]] <- abs(sqrt(summary(list.reg[[i]])$"resVar") / mean(fitted(list.reg[[i]])))
    # Plot
    plot(list.reg[[i]], col=col.wheel[i], 
         xlab="MAb Concentration (ng/ml)",
         ylab="% Infection", ylim=c(0,1.2), 
         pch=19, cex=0.5)
    arrows(conc, cor.st[1,], conc, cor.st[2,], length=0.05, angle=90, code=3, lty=3, col=col.wheel[i])
    par(new=T)
  }}
legend("right", "top", legend=abs.left, fill=col.wheel[1:length(abs.right)]) 
# Name the listss 
names(list.reg) <- abs.left
names(list.eds) <- abs.left
# Plot the table of EDs and %GFP at max antibody for each 
mat.plot <- as.data.frame(matrix(nrow=length(list.reg), ncol=5))
colnames(mat.plot) <- c("ab ID", "ED50", "ED50 error", "%GFP at ab max", "RSD")
for(i in 1:length(list.reg)){
  mat.plot[i,1] <- abs.left[i]
  mat.plot[i,2] <- round(list.eds[[i]][3,1],2) # get ED 50 
  mat.plot[i,3] <- round(list.eds[[i]][3,2],2) # get ED 50 st. error 
  mat.plot[i,4] <- round(list.gfp.max[[i]],2)
  mat.plot[i,5] <- round(list.rsd[[i]],3)
}
# Plot table 
tbl <- tableGrob(mat.plot)
plot(tbl)
dev.off()
# # 
message("Making right plate pdf")
pdf(paste(args[3], args[2], "_right_plate_reg.pdf", sep=""),10,10)
# Iterate for every row pair 
for(i in 1:nrow(pairs)){
  # Subset data from main data frame 
  pair.1 <- data[grep(pairs[i,1], data$Well.ID),][11:20,] # 11:20 selects X14-X23
  pair.2 <- data[grep(pairs[i,2], data$Well.ID),][11:20,]
  # Get GFP value and normalize to no antibody control (last entry X23)
  # grep for 23 to get the X23 column from both replicates
  con.pair.1 <- pair.1$X.GFP.[grep(23, pair.1$Well.ID)] 
  con.pair.2 <- pair.2$X.GFP.[grep(23, pair.2$Well.ID)]
  # Correct the GFP values by dividing with the no antibody control 
  cor.pair.1.t <- pair.1$X.GFP./con.pair.1
  cor.pair.2.t <- pair.2$X.GFP./con.pair.2
  # Remove the zero value 
  cor.pair.1 <- cor.pair.1.t[-length(cor.pair.1.t)]
  cor.pair.2 <- cor.pair.2.t[-length(cor.pair.2.t)]
  # Get average values 
  cor.means <- rowMeans(cbind(cor.pair.1, cor.pair.2))
  cor.st <- apply(cbind(cor.pair.1, cor.pair.2), MARGIN = 1, range)
  # DOES IT MAKE SENSE ? 
  if(length(which(is.infinite(cor.means)==T)) ==0){
    # Fit using a 4-parameter log logistic model and store the model
    list.reg[[i]] <- drm(cor.means ~ conc,fct = LL.4())
    # Get ED5,10 and 50 values 
    list.eds[[i]] <- ED(object = list.reg[[i]], c(5,10,50), interval = "delta", display=F)
    # Also store %GFP and max antibody 
    list.gfp.max[[i]] <- cor.means[1]*100
    # Also store the residual standard error 
    list.rsd[[i]] <- abs(sqrt(summary(list.reg[[i]])$"resVar") / mean(fitted(list.reg[[i]])))
    # Plot
    plot(list.reg[[i]], col=col.wheel[i], 
         xlab="MAb Concentration (ng/ml)",
         ylab="% Infection", ylim=c(0,1.2), 
         pch=19, cex=0.5)
    arrows(conc, cor.st[1,], conc, cor.st[2,], length=0.05, angle=90, code=3, lty=3, col=col.wheel[i])
    par(new=T)
  }}
legend("right", "top", legend=abs.right, fill=col.wheel[1:length(abs.right)]) 
# Name the listss 
names(list.reg) <- abs.right
names(list.eds) <- abs.right
# Plot the table of EDs and %GFP at max antibody for each 
mat.plot <- as.data.frame(matrix(nrow=length(list.reg), ncol=5))
colnames(mat.plot) <- c("ab ID", "ED50", "ED50 error", "%GFP at ab max", "RSD")
for(i in 1:length(list.reg)){
  mat.plot[i,1] <- abs.right[i]
  mat.plot[i,2] <- round(list.eds[[i]][3,1],2) # get ED 50 
  mat.plot[i,3] <- round(list.eds[[i]][3,2],2) # get ED 50 st. error 
  mat.plot[i,4] <- round(list.gfp.max[[i]],2)
  mat.plot[i,5] <- round(list.rsd[[i]],3)
}
# Plot table 
tbl <- tableGrob(mat.plot)
plot(tbl)
dev.off()