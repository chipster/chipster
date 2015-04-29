# TOOL create_cloud_project.R: "Create new cloud session." (Create new cloud session from server side fastq and bam files. An error box will pop up when this tool has finished. Open the Details and read the first line for further information.)
# PARAMETER project.name: "The project name" TYPE STRING (Enter a name for the new cloud project)
# PARAMETER OPTIONAL file.filter: "A filter string" TYPE STRING (only files containing this string are imported: *file.filter*)
# PARAMETER OPTIONAL file.extension: "The file extension" TYPE STRING (only files ending with file.extension are imported : *file.extension)

# OH 29.04.2015

user<-chipster.R.job.calling.user
project<-project.name
filter<-file.filter
extension<-file.extension

filter<-paste("filter=",filter,sep<-"",collapse<-"")
filter<-gsub("\\s+","",filter)
extension<-paste("extension=",extension,sep<-"",collapse<-"")
extension<-gsub("\\s+","",extension)

binary<-file.path(chipster.module.path ,"/shell/create_cloud_project.sh")
binary.full<-paste(binary,user,project,filter,extension," > create_cloud_project.log 2>&1",sep<-" ")

if( nchar(user) <= 0 || nchar(project) <= 0 || system(binary.full) != 0 ) {
	out<-"The cloud project can not be created\r\n\r\n"
} else {
	out<-"The cloud project has been successful created. You can open it using File -> Open cloud session...\r\n\r\n"
}

log<-paste(readLines("create_cloud_project.log"),sep="",collapse="\r\n")

out<-paste(out,"log file:\r\n\r\n",log,sep="")

stop(out)

