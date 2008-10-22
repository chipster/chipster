Quick installation notes for running Chipster service in a single Linux system.


Requirements
============
- Java 1.5+
- R-2.6.1
- the following tcp ports open in the firewall
	- 61616 for message broker service
	- 8080 for file broker service
	- 8081 for webstart service (optional)


1) Installing needed R libraries
================================
To install the R libraries needed by Chipster, run (as root if needed):

	./install_r_pkg.sh


2) Configuring Chipster services
================================
To configure the Chipster services, run the following to scripts. Both 
scripts will ask for confirmation before writing changes to files. 
Defaults should be fine for a local installation. 

	./configure.sh
	./genpasswd.sh
	

3) Starting and stopping Chipster services
==========================================
To start all the Chipster services, run:

	./chipster start
	
In addition to 'start', 'stop', 'restart', and 'status' can also be used.


4) Testing the installation with the client
===========================================
To start the client locally (on the same machine as the services), run:

	./client/bin/chipster-client

To start the client from other machines using Java Web Start, go to the 
Web Start address specified when running the configure.sh. Default address is:

	http://hostname-of-this-machine:8081 

The default username/password is chipster/chipster.
	