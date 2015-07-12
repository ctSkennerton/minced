JFLAGS = -g
JC = javac
JAR = jar
.SUFFIXES: .java .class

.java.class:
	$(JC) $(JFLAGS) $*.java

CLASSES = CRISPR.java CRISPRFinder.java CRISPRUtil.java DNASequence.java FASTAReader.java SearchUtil.java minced.java IntervalSearchTree.java

default: classes minced.jar

classes: $(CLASSES:.java=.class)

minced.jar: classes
	$(JAR) cfm minced.jar MANIFEST.txt *class

clean:
	$(RM) *.class

test: minced.jar
	@echo "Testing..."
	@./minced -gff t/Aquifex_aeolicus_VF5.fna >t/Aquifex_aeolicus_VF5.output
	@diff -q t/Aquifex_aeolicus_VF5.output t/Aquifex_aeolicus_VF5.expected > /dev/null || (echo "Failed" && exit 1)
	@echo "Passed"

