# -------- Compile Socrates ---------
socrates=`dirname $0`/
libs=${socrates}lib/sam-1.77.jar:${socrates}lib/commons-lang3-3.1.jar:${socrates}lib/picard-1.85.jar:${socrates}lib/snappy-java-1.0.3-rc3.jar:${socrates}lib/commons-cli-1.2.jar

echo Compiling classes
javac -d ${socrates}bin -sourcepath ${socrates}src -cp $libs ${socrates}src/net/wehi/socrates/*.java ${socrates}src/net/wehi/socrates/util/*.java
if [[ $? == 1 ]]
then
    echo "Compile failed."
    exit 1
fi

