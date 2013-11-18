JFLAGS = -g
JC = javac
JAR = jar
.SUFFIXES: .java .class

.java.class:
	$(JC) $(JFLAGS) $*.java

CLASSES = CRISPR.java CRISPRFinder.java CRISPRUtil.java DNASequence.java FASTAReader.java SearchUtil.java minced.java

default: classes minced

classes: $(CLASSES:.java=.class)

minced: classes
	$(JAR) cfm minced.jar MANIFEST.txt *class

clean:
	$(RM) *.class
