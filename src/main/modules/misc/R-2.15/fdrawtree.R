# TOOL fdrawtree.R: "Plot an unrooted tree" (Plots an unrooted tree diagram)
# INPUT OPTIONAL intreefile: intreefile TYPE GENERIC 
# OUTPUT OPTIONAL tree.ps
# OUTPUT OPTIONAL tree.pdf
# OUTPUT OPTIONAL fdrawtree.log
# PARAMETER OPTIONAL iterate: "Iterate to improve tree" TYPE [n: No, e: "Equal-Daylight algorithm", b: "n-Body algorithm"] FROM 1 TO 1 DEFAULT e (Iterate to improve tree)
# PARAMETER OPTIONAL lengths: "Use branch lengths from user trees" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Use branch lengths from user trees)
# PARAMETER OPTIONAL labeldirection: "Label direction" TYPE [a: along, f: fixed, r: radial, m: middle] FROM 1 TO 1 DEFAULT m (Label direction)
# PARAMETER OPTIONAL treeangle: "Angle the tree is to be plotted" TYPE DECIMAL FROM -360.0 TO 360.0 DEFAULT 90.0 (Angle the tree is to be plotted)
# PARAMETER OPTIONAL arc: "Degrees the arc should occupy" TYPE DECIMAL FROM 0.0 TO 360.0 DEFAULT 360 (Degrees the arc should occupy)
# PARAMETER OPTIONAL labelrotation: "Angle of labels" TYPE DECIMAL FROM 0.0 TO 360.0 DEFAULT 90.0 (Angle of labels (0 degrees is horizontal for a tree growing vertically\))
# PARAMETER OPTIONAL rescaled: "Automatically rescale branch lengths" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Automatically rescale branch lengths)
# PARAMETER OPTIONAL bscale: "Centimeters per unit branch length" TYPE DECIMAL DEFAULT 1.0 (Centimeters per unit branch length)
# PARAMETER OPTIONAL treedepth: "Depth of tree as fraction of its breadth" TYPE DECIMAL FROM 0.1 TO 100.0 DEFAULT 0.53 (Depth of tree as fraction of its breadth)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")


emboss.binary <- file.path(emboss.path, "fdrawtree")
emboss.parameters <- paste("-auto -plotter l -previewer n")
emboss.parameters <- paste(emboss.parameters, "-intreefile intreefile")
emboss.parameters <- paste(emboss.parameters, "-plotfile tree.ps")
emboss.parameters <- paste(emboss.parameters, "-iterate", iterate)
emboss.parameters <- paste(emboss.parameters, "-lengths", lengths)
emboss.parameters <- paste(emboss.parameters, "-labeldirection", labeldirection)
emboss.parameters <- paste(emboss.parameters, "-treeangle", treeangle)
emboss.parameters <- paste(emboss.parameters, "-arc", arc)
emboss.parameters <- paste(emboss.parameters, "-labelrotation", labelrotation)
emboss.parameters <- paste(emboss.parameters, "-rescaled", rescaled)
emboss.parameters <- paste(emboss.parameters, "-bscale", bscale)
emboss.parameters <- paste(emboss.parameters, "-treedepth", treedepth)


command.full <- paste(emboss.binary, emboss.parameters, ' >> fdrawtree.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> fdrawtree.log' )
system(echo.command)

system(command.full)

system("ps2pdf tree.ps")
system("ls -l >> fdrawtree .log")

if ( save_log == "no") {
	system ("rm -f fdrawtree.log")
}
