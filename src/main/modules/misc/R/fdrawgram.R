# TOOL fdrawgram.R: "Plot a cladogram" (Plots a cladogram- or phenogram-like rooted tree diagram using PHYLIP program drawgam)
# INPUT OPTIONAL intreefile: intreefile TYPE GENERIC 
# OUTPUT OPTIONAL cladocram.ps
# OUTPUT OPTIONAL cladocram.pdf
# OUTPUT OPTIONAL fdrawgram.log
# PARAMETER OPTIONAL fontfile: "Fontfile name" TYPE STRING DEFAULT font1 (Fontfile name)
# PARAMETER OPTIONAL grows: "Tree grows horizontally" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Tree grows horizontally)
# PARAMETER OPTIONAL style: "Tree style output" TYPE [c: "cladogram (v-shaped\)", p: "phenogram (branches are square\)", v: "curvogram (branches are 1/4 out of an ellipse\)", e: "eurogram (branches angle outward, then up\)", s: "swooporam (branches curve outward then reverse\)", o: "circular tree"] FROM 1 TO 1 DEFAULT c (Tree style output)
# PARAMETER OPTIONAL lengths: "Use branch lengths from user trees" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Use branch lengths from user trees)
# PARAMETER OPTIONAL labelrotation: "Angle of labels (0 degrees is horizontal for a tree growing vertically\)" TYPE DECIMAL FROM 0.0 TO 360.0 DEFAULT 90.0 (Angle of labels (0 degrees is horizontal for a tree growing vertically\))
# PARAMETER OPTIONAL rescaled: "Automatically rescale branch lengths" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Automatically rescale branch lengths)
# PARAMETER OPTIONAL bscale: "Centimeters per unit branch length" TYPE DECIMAL DEFAULT 1.0 (Centimeters per unit branch length)
# PARAMETER OPTIONAL treedepth: "Depth of tree as fraction of its breadth" TYPE DECIMAL FROM 0.1 TO 100.0 DEFAULT 0.53 (Depth of tree as fraction of its breadth)
# PARAMETER OPTIONAL stemlength: "Stem length as fraction of tree depth" TYPE DECIMAL FROM 0.01 TO 100.0 DEFAULT 0.05 (Stem length as fraction of tree depth)
# PARAMETER OPTIONAL nodespace: "Character height as fraction of tip spacing" TYPE DECIMAL FROM 0.1 TO 100.0 DEFAULT 0.3333 (Character height as fraction of tip spacing)
# PARAMETER OPTIONAL nodeposition: "Position of interior nodes" TYPE [i: "Intermediate between their immediate descendants", w: "Weighted average of tip positions", c: "Centered among their ultimate descendants", n: "Innermost of immediate descendants", v: "So tree is v shaped"] FROM 1 TO 1 DEFAULT c (Position of interior nodes)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

options(scipen=999)
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

emboss.binary <- file.path(emboss.path, "fdrawgram")
emboss.parameters <- paste("-auto -plotter l -previewer n")
emboss.parameters <- paste(emboss.parameters, "-intreefile intreefile" )
emboss.parameters <- paste(emboss.parameters, "-plotfile cladocram.ps")
emboss.parameters <- paste(emboss.parameters, "-grows", grows)
emboss.parameters <- paste(emboss.parameters, "-style", style)
emboss.parameters <- paste(emboss.parameters, "-lengths", lengths)
emboss.parameters <- paste(emboss.parameters, "-labelrotation", labelrotation)
emboss.parameters <- paste(emboss.parameters, "-rescaled", rescaled)
emboss.parameters <- paste(emboss.parameters, "-bscale", bscale)
emboss.parameters <- paste(emboss.parameters, "-treedepth", treedepth)
emboss.parameters <- paste(emboss.parameters, "-stemlength", stemlength)
emboss.parameters <- paste(emboss.parameters, "-nodespace", nodespace)
emboss.parameters <- paste(emboss.parameters, "-nodeposition", nodeposition)


command.full <- paste(emboss.binary, emboss.parameters, ' >> fdrawgram.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> fdrawgram.log' )
system(echo.command)

#stop(paste('CHIPSTER-NOTE:  ', command.full ))

system(command.full)

system("ps2pdf cladocram.ps")
system("ls -l >> fdrawgram.log")

if ( save_log == "no") {
	system ("rm -f fdrawgram.log")
}
