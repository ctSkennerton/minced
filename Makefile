JFLAGS = -g 
JC = javac
JAR = jar
.SUFFIXES: .java .class
JAVA_FILES = CRISPR.java CRISPRFinder.java CRISPRUtil.java DNASequence.java FASTAReader.java SearchUtil.java minced.java IntervalSearchTree.java
JAVA_CLASSES = $(JAVA_FILES:.java=.class)
.PHONY: default clean test

.java.class:
	$(JC) $(JFLAGS) $*.java

default: minced.jar

minced.jar: $(JAVA_CLASSES)
	$(JAR) cfm minced.jar MANIFEST.txt *.class

clean:
	$(RM) $(JAVA_CLASSES)

test: minced.jar
	@echo "Testing..."
	@./minced -gff t/Aquifex_aeolicus_VF5.fna >t/Aquifex_aeolicus_VF5.output
	@diff -q t/Aquifex_aeolicus_VF5.output t/Aquifex_aeolicus_VF5.expected > /dev/null || (echo "Failed" && exit 1)
	@echo "Passed"

