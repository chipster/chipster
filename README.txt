Chipster Source Package
-----------------------

Out-of-the-box this project can be built with Ant and used with Eclipse.
Source package structure follows Maven layout, but this is not a Maven project. 
Unlike in Maven, xternal dependencies are located in "ext" directory.

For license, see LICENSE.TXT.


Building
--------

For building options, use "ant -p" in project root directory.

For compiling all Java sources, use "ant compile".

For complete build, use "ant". Please note that you need to have keystore 
available because client JAR needs to be signed. Any key will do, so 
you can use "keytool" to generate your signing key.

To build a freshly checked out project directory, go to command line and
in the project directory give following command:

keytool -genkey -alias csc -keystore keystore.ks -storepass chipster -validity 1825

Fill in your data to keytool. You might want to use your organisations name 
also as your own name as it is shown by Java Web Start. Use the keystore password
also as a key password.

Next command "ant", and the project will be built and packages are available in "dist"
directory.



