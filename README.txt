Chipster Source Package
-----------------------

For compiling all Java sources, use "ant compile".

For complete build, use "ant". Please note that you need to have keystore 
available because client JAR needs to be signed. Any key will do, so 
you can use "keytool -genkey" to generate your signing key.

For building options, use "ant -p" in package root directory.

Out-of-the-box this project can be built with Ant and used with Eclipse.
Source package structure follows Maven layout, but this is not a Maven project. 
Especially external dependencies are located in "ext" directory.

For license, see LICENSE.TXT.
