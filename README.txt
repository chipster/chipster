Chipster v3 has reached end-of-life
-----------------------

Up to version 3.x Chipster was a desktop application requiring Java Web Start. 
Chipster version 4 is a Web application which runs on your browser. 

Here are some links to get started with the new Chipster v4:
- Main website: [https://chipster.csc.fi]
- Installation instructions: [https://github.com/chipster/chipster-openshift/tree/k3s/k3s]
- Source code: [https://github.com/chipster]

Chipster Source Package
-----------------------

Out-of-the-box this project can be built with Ant and used with Eclipse.
Source package structure follows Maven layout, but this is not a Maven project. 
Unlike in Maven, xternal dependencies are located in "ext" directory.

This is a version of Chipster which is supposed to be running with Chipster Jog Manager component.

For license, see LICENSE.TXT.


Building
--------

For building options, use "ant -p" in project root directory.

For complete build, use "ant". Please note that you need to have keystore 
available because client JAR needs to be signed. Any key will do, so 
you can use "keytool" to generate your signing key.

To build a freshly checked out project directory, go to command line and
in the project directory give following command:

keytool -genkey -alias chipster -keystore keystore.ks -storepass chipster -validity 1825

Fill in your data to keytool command parameters. You might want to use your organisations name 
also as your own name as it is shown by Java Web Start. Use the keystore password
also as a key password.

Next, issue command "ant", and the project will be built and packages will be available in "dist"
directory.

If you do not want to place your keystore in the project directory, you can create
a file called alternative-keystore-path.txt and place path to your keystore into it.
Do not enter a newline after the path! The build script will read that file and
use the alternative path into keystore. This is strongly recommended if you do also 
commits into the public version repository, because placing a private keystore file into 
project directory creates a risk of publishing it by accident.
