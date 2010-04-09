# ANALYSIS Utilities/"Intersect lists" (Performs Boolean operations on 2 or 3 data tables or lists
# - that have one column in common - to identify the rows that intersect, are unique or form a union. The
# results are collected in one single table with columns for each of the operations. A Venn diagram giving
# a visual interpretation of the results is also returned.)
# INPUT GENERIC genelist[...].tsv OUTPUT intersect-lists-operation.tsv, intersect-lists-summary.tsv, venn-diagram-plot.png
# PARAMETER common.column STRING DEFAULT empty (The name of the column that is common to the data tables.)
# PARAMETER operation [intersection, union] DEFAULT intersection (Defines the operation to be performed on the
# lists, where "intersection" will yield a list of only the rows that are common between lnput lists, whereas "union"
# will yield a list of the rows that appear in any of the input lists.)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# INPUT GENERIC genelist[...].tsv OUTPUT intersect-lists-operation.tsv, intersect-lists-summary.tsv, venn-diagram-plot.png
#common.column <- "symbol"
#operation <- "intersect"
#image.width <- 5
#image.height <- 5
h <- image.height
w <- image.width 

# Sanity check the input data
files<-dir()
files<-files[grep("genelist", files)]
number_files <- length(files)
if (number_files > 3 | number_files < 2) {
	stop("CHIPSTER-NOTE: You need to have 2 or 3 input files to use this tool!")
}

# Read in data files
for (count in 1:number_files) {
	assign (paste("data_", count, sep=""), read.table(files[count], header=T, sep="\t")) 
}

# Extract the common columns from the input data and condense to unique values
# to get rid of multiple entries
for (count in 1:number_files) {
	data_temp <- get(paste("data_", count, sep=""))
	assign (paste("list_", count, sep=""), unique (data_temp [,names(data_temp)==common.column]))
}

# remove "NA" entries
if (number_files==2) {
	list_1 <- list_1[-is.na(list_1)]
	list_2 <- list_2[-is.na(list_2)]
}
if (number_files==3) {
	list_1 <- list_1[-is.na(list_1)]
	list_2 <- list_2[-is.na(list_2)]
	list_3 <- list_3[-is.na(list_3)]	
}

if (number_files==2) {
	
	# find out intersection between all 3 lists
	intersect_1_2 <- intersect(list_1, list_2)
	if (operation=="intersection") {
		result_operation <- intersect_1_2
	}	
	
	# find out union between pairs
	union_1_2 <- union(list_1, list_2)
	if (operation=="union") {
		result_operation <- union_1_2
	}	
	
	# find out unique for each list
	unique_1 <- setdiff (list_1, intersect_1_2)
	unique_2 <- setdiff (list_2, intersect_1_2)
	
	# set up plotting area
	bitmap(file="venn-diagram-plot.png", width=w/72, height=h/72)
#	bmp(file="venn-diagram-plot.png", width=w, height=h, res=72, units="px")
#	png(file="venn-diagram-plot.png")
	plot(-1:1, -1:1, type="n", axes = FALSE, xlab = "", ylab = "")
	
	# draw overlapping circles
	symbols(-.25,0, circles = .5, inches = FALSE, add = T, fg=3)
	symbols(.25,0, circles = .5, inches = FALSE, add = T, fg=4)
	
	# place list identifiers including length numbers 
	text(-.4,0.75, labels=paste("list_1","\n",length(list_1)), col="black", cex=1.2)
	text(.4,0.75, labels=paste("list_2","\n",length(list_2)), col="black", cex=1.2)
	
	# place length numbers for the intersections
	text(0,0, labels=length(intersect_1_2), col="black", cex=1.2)
	
	# place length numbers for the unique non-overlapping areas
	text(-.4,0, labels=length(unique_1), col="black", cex=1.2)
	text(.4,0, labels=length(unique_2), col="black", cex=1.2)
	
	# shut down plotting device
	dev.off()
	
	# output the result table
	longest_list <- length(union_1_2)
	result_table <- matrix(nrow=longest_list,ncol=6)
	list_1 <- append(as.character(list_1),rep("",longest_list-length(list_1)))
	list_2 <- append(as.character(list_2),rep("",longest_list-length(list_2)))
	unique_1 <- append(as.character(unique_1),rep("",longest_list-length(unique_1)))
	unique_2 <- append(as.character(unique_2),rep("",longest_list-length(unique_2)))
	intersect_1_2 <- append(as.character(intersect_1_2),rep("",longest_list-length(intersect_1_2)))
	union_1_2 <- append(as.character(union_1_2),rep("",longest_list-length(union_1_2)))
	result_table [,1] <- as.character(list_1)
	result_table [,2] <- as.character(list_2)
	result_table [,3] <- as.character(unique_1)
	result_table [,4] <- as.character(unique_2)
	result_table [,5] <- as.character(intersect_1_2)
	result_table [,6] <- as.character(union_1_2)
	rownames(result_table) <- as.character(1:longest_list)
	colnames(result_table) <- c(
			"list_1",
			"list_2",
			"unique_1",
			"unique_2",
			"intersect_1_2",
			"union_1_2")
	result_operation <- as.data.frame(result_operation)
	colnames(result_operation) <- common.column
	write.table(result_operation, file="intersect-lists-operation.tsv", sep="\t", row.names=T, col.names=T, quote=F)
	write.table(result_table, file="intersect-lists-summary.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

if (number_files==3) {
	
	# find out intersection between all 3 lists
	intersect_1_2_3 <- intersect(intersect(list_1, list_2), list_3)
	if (operation=="intersection") {
		result_operation <- intersect_1_2_3
	}	
	
	# find out intersections between pairs
	intersect_1_2 <- intersect(list_1, list_2)
	intersect_1_3 <- intersect(list_1, list_3)
	intersect_2_3 <- intersect(list_2, list_3)
	intersect_1_2_plot <- setdiff(intersect(list_1, list_2), intersect_1_2_3)
	intersect_1_3_plot <- setdiff(intersect(list_1, list_3), intersect_1_2_3)
	intersect_2_3_plot <- setdiff(intersect(list_2, list_3), intersect_1_2_3)
	
	# find out union between pairs
	union_1_2 <- union(list_1, list_2)
	union_1_3 <- union(list_1, list_3)
	union_2_3 <- union(list_2, list_3)
	
	# find out the union between all 3 lists
	union_1_2_3 <- union(union(union(union_1_2, union_1_3), union_2_3), intersect_1_2_3)
	if (operation=="union") {
		result_operation <- union_1_2_3
	}	
	
	# find out unique for each list
	unique_1 <- setdiff (list_1, union(union(intersect_1_2, intersect_1_3), intersect_1_2_3))
	unique_2 <- setdiff (list_2, union(union(intersect_1_2, intersect_2_3), intersect_1_2_3))
	unique_3 <- setdiff (list_3, union(union(intersect_1_3, intersect_2_3), intersect_1_2_3))
	
	# set up plotting area
bitmap(file="venn-diagram-plot.png", width=w/72, height=h/72)
#	bmp(file="venn-diagram-plot.png", width=w, height=h, res=72, units="px")
#	bmp(file="venn-diagram-plot.png")
	plot(-1:1, -1.3:1, type="n", axes = FALSE, xlab = "", ylab = "")
	
	# draw overlapping circles
	symbols(0,-.5, circles = .5, inches = FALSE, add = T, fg=2)
	symbols(-.25,0, circles = .5, inches = FALSE, add = T, fg=3)
	symbols(.25,0, circles = .5, inches = FALSE, add = T, fg=4)
	
	# place list identifiers including length numbers 
	text(-.9,.6, labels=paste("list_1","\n",length(list_1)), col="black", cex=1.2)
	text(.8,.6, labels=paste("list_2","\n",length(list_2)), col="black", cex=1.2)
	text(0,-1.25, labels=paste("list_3","\n",length(list_3)), col="black", cex=1.2)
	
	# place length numbers for the intersections
	text(0,-.15, labels=length(intersect_1_2_3), col="black", cex=1.2)
	text(-.3,-.35, labels=length(intersect_1_3_plot), col="black", cex=1.2)
	text(.3,-.35, labels=length(intersect_2_3_plot), col="black", cex=1.2)
	text(0,.2, labels=length(intersect_1_2_plot), col="black", cex=1.2)
	
	# place length numbers for the unique non-overlapping areas
	text(0,-.8, labels=length(unique_3), col="black", cex=1.2)
	text(-.4,.3, labels=length(unique_1), col="black", cex=1.2)
	text(.4,.3, labels=length(unique_2), col="black", cex=1.2)
	
	# place length numbers for the pairwise unions 
	text(0,.7, labels=paste("Union_1_2","\n",length(union_1_2)), col="black", cex=1.2)
	text(-.85,-.65, labels=paste("Union_1_3","\n",length(union_1_3)), col="black", cex=1.2)
	text(.85,-.65, labels=paste("Union_2_3","\n",length(union_2_3)), col="black", cex=1.2)
	
	# shut down plotting device
	dev.off()
	
	# output the result table
	longest_list <- length(union_1_2_3)
	result_table <- matrix(nrow=longest_list,ncol=14)
	list_1 <- append(as.character(list_1),rep("",longest_list-length(list_1)))
	list_2 <- append(as.character(list_2),rep("",longest_list-length(list_2)))
	list_3 <- append(as.character(list_3),rep("",longest_list-length(list_3)))
	unique_1 <- append(as.character(unique_1),rep("",longest_list-length(unique_1)))
	unique_2 <- append(as.character(unique_2),rep("",longest_list-length(unique_2)))
	unique_3 <- append(as.character(unique_3),rep("",longest_list-length(unique_3)))
	intersect_1_2 <- append(as.character(intersect_1_2),rep("",longest_list-length(intersect_1_2)))
	intersect_1_3 <- append(as.character(intersect_1_3),rep("",longest_list-length(intersect_1_3)))
	intersect_2_3 <- append(as.character(intersect_2_3),rep("",longest_list-length(intersect_2_3)))
	intersect_1_2_3 <- append(as.character(intersect_1_2_3),rep("",longest_list-length(intersect_1_2_3)))
	union_1_2 <- append(as.character(union_1_2),rep("",longest_list-length(union_1_2)))
	union_1_3 <- append(as.character(union_1_3),rep("",longest_list-length(union_1_3)))
	union_2_3 <- append(as.character(union_2_3),rep("",longest_list-length(union_2_3)))
	result_table [,1] <- as.character(list_1)
	result_table [,2] <- as.character(list_2)
	result_table [,3] <- as.character(list_3)
	result_table [,4] <- as.character(unique_1)
	result_table [,5] <- as.character(unique_2)
	result_table [,6] <- as.character(unique_3)
	result_table [,7] <- as.character(intersect_1_2)
	result_table [,8] <- as.character(intersect_1_3)
	result_table [,9] <- as.character(intersect_2_3)
	result_table [,10] <- as.character(intersect_1_2_3)
	result_table [,11] <- as.character(union_1_2)
	result_table [,12] <- as.character(union_1_3)
	result_table [,13] <- as.character(union_2_3)
	result_table [,14] <- as.character(union_1_2_3)
	rownames(result_table) <- as.character(1:longest_list)
	colnames(result_table) <- c(
			"list_1",
			"list_2",
			"list_3",
			"unique_1",
			"unique_2",
			"unique_3",
			"intersect_1_2",
			"intersect_1_3",
			"intersect_2_3",
			"intersect_all",		
			"union_1_2",
			"union_1_3",
			"union_2_3",
			"union_all")
	result_operation <- as.data.frame(result_operation)
	colnames(result_operation) <- common.column
	write.table(result_operation, file="intersect-lists-operation.tsv", sep="\t", row.names=T, col.names=T, quote=F)
	write.table(result_table, file="intersect-lists-summary.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
