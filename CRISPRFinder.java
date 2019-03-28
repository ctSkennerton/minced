import java.io.*;
import java.util.Vector;

public class CRISPRFinder
{
   private String inputFileName;
   private String outputFileName;
   private String outputGffFileName;

   private int screenDisplay;
   private int minNumRepeats;
   private int minRepeatLength;
   private int maxRepeatLength;
   private int minSpacerLength;
   private int maxSpacerLength;
   private int searchWindowLength;
   private int outputformat;
   private boolean printSpacers;
   private PrintStream spacers;
   private boolean printGffHeader;

   DNASequence sequence = null;
   int sequenceLength = 0;
   private int totalCrisprCount = 0;


   public CRISPRFinder(String _inputFileName,
		   String _outputFileName,
		   String _outputGffFileName,
		   int _screenDisplay,
		   int _minNumRepeats,
		   int _minRepeatLength,
		   int _maxRepeatLength,
		   int _minSpacerLength,
		   int _maxSpacerLength,
		   int _searchWindowLength,
		   int _outputformat,
		   boolean _spacers)
   {
      inputFileName = _inputFileName;
      outputFileName = _outputFileName;
      outputGffFileName = _outputGffFileName;

      screenDisplay = _screenDisplay;
      minNumRepeats = _minNumRepeats;
      minRepeatLength = _minRepeatLength;
      maxRepeatLength = _maxRepeatLength;
      minSpacerLength = _minSpacerLength;
      maxSpacerLength = _maxSpacerLength;
      searchWindowLength = _searchWindowLength;
      outputformat = _outputformat;
      printSpacers = _spacers;
      printGffHeader = true;

      if(_spacers) {
    	  File outputFile;
    	  FileOutputStream outputFileStream;

    	  if(_outputFileName == "") {
        	  String input_file_prefix = removeExtention(inputFileName);
    		  outputFile = new File(input_file_prefix + "_spacers.fa");

    	  } else {
        	  String input_file_prefix = removeExtention(outputFileName);
        	  outputFile = new File(input_file_prefix + "_spacers.fa");
    	  }

    	  try {
    		    outputFileStream = new FileOutputStream(outputFile, false);
        	  spacers = new PrintStream(outputFileStream);
    	  } catch (FileNotFoundException e) {
    		  // TODO Auto-generated catch block
    		  e.printStackTrace();
    	  }
      }

   }

   private static String removeExtention(String filePath) {
	    // These first few lines the same as Justin's
	    File f = new File(filePath);

	    // if it's a directory, don't remove the extention
	    if (f.isDirectory()) return filePath;

	    String name = f.getName();

	    // Now we know it's a file - don't need to do any special hidden
	    // checking or contains() checking because of:
	    final int lastPeriodPos = name.lastIndexOf('.');
	    if (lastPeriodPos <= 0)
	    {
	        // No period after first character - return name as it was passed in
	        return filePath;
	    }
	    else
	    {
	        // Remove the last period and everything after it
	        File renamed = new File(f.getParent(), name.substring(0, lastPeriodPos));
	        return renamed.getPath();
	    }
	}

   public boolean goCRISPRFinder()
   {
      File inputFile = new File(inputFileName);
      if (! inputFile.exists())
      {
         System.err.println("File name does not exist:  " + inputFile.getPath());
         return false;
      }
      //System.out.println("\n\nReading file '" + inputFile.getPath() + "'");
      FASTAReader fastaReader = new FASTAReader(inputFile.getPath());
      if (! fastaReader.isFASTA())
      {
         System.err.println("File '" + inputFile.getPath() + "' does not seem to be FASTA-formatted.");
         return false;
      }

      while ( fastaReader.read() )
      {
         sequence = new DNASequence(fastaReader.getSequence(), fastaReader.getName(), fastaReader.getDesc());
         //System.out.println("Sequence '" + sequence.getName() + "' (" + sequence.length() + " bp)");

         /*
         // Check that sequence is a valid DNA sequence
         if ( ! sequence.isDNASequence() )
         {
            System.out.println(sequence.getErrorLog() + "Skipping this sequence!");
            System.out.println("");
            continue;
         }
         */
         sequence.mask(100);

         try
         {
            findRepeats(sequence, fastaReader.getNum());
         }
         catch (Exception e)
         {
            System.err.println ("Error processing input file '" + inputFile.getPath() + "'. Please, check contents.\n");
            e.printStackTrace(System.err);
         }

      }

      //System.out.println("Done!");
      return true;
   }


   private boolean findRepeats( DNASequence sequence, int readNum )
   {

      Vector<CRISPR> CRISPRVector = new Vector<CRISPR>();
      sequenceLength = sequence.length();
      int actualRepeatLength;
      boolean repeatsFound = false;

      CRISPR candidateCRISPR;
      String pattern;

      if ((searchWindowLength < 6) || (searchWindowLength > 9))
      {
         // Change window length
         int oldSearchWindowlength = searchWindowLength;
         searchWindowLength = 8;
         System.err.println("Changing window length to " + searchWindowLength + " instead of " + oldSearchWindowlength);
      }

      double spacerToSpacerMaxSimilarity = 0.62;
      int spacerToSpacerLengthDiff = 12;
      int spacerToRepeatLengthDiff = 30;

      //the mumber of bases that can be skipped while we still guarantee that the entire search
      //window will at some point in its iteration thru the sequence will not miss a any repeat
      int skips = minRepeatLength - (2 * searchWindowLength - 1);
      if (skips < 1)
         skips = 1;

      //System.out.println("Searching for repeats...");
      long repeatSearchStart = System.currentTimeMillis();

      SearchUtil searchUtil = new SearchUtil();

      int searchEnd = sequenceLength - maxRepeatLength - maxSpacerLength - searchWindowLength;
      for (int j = 0; j <= searchEnd; j = j + skips)
      {
         candidateCRISPR = new CRISPR();

         int beginSearch = j + minRepeatLength + minSpacerLength;
         int endSearch = j + maxRepeatLength + maxSpacerLength + searchWindowLength;

         if (endSearch > sequenceLength)
            endSearch = sequenceLength;

         if (endSearch < beginSearch)  //should never occur
            endSearch = beginSearch;

        int mask_end = sequence.masked(beginSearch, endSearch);
        if(mask_end != -1)
        {
            j = mask_end;
            continue;
        }

         String text = sequence.substring(beginSearch, endSearch);
         pattern = sequence.substring(j, j + searchWindowLength);

         //if pattern is found, add it to candidate list and scan right for additional similarly spaced repeats
         int patternInTextIndex = searchUtil.boyer_mooreSearch(text, pattern);
         if (patternInTextIndex >= 0)
         {
            candidateCRISPR.addRepeat(j);
            candidateCRISPR.addRepeat(beginSearch + patternInTextIndex);
            candidateCRISPR = CRISPRUtil.scanRight(candidateCRISPR, pattern, minSpacerLength, 24, searchUtil);
         }

         if ( (candidateCRISPR.numRepeats() >= minNumRepeats) )  //make sure minNumRepeats is always at least 2
         {
            candidateCRISPR = CRISPRUtil.getActualRepeatLength(candidateCRISPR, searchWindowLength, minSpacerLength, minNumRepeats);
            actualRepeatLength = candidateCRISPR.repeatLength();

            if ( (actualRepeatLength >= minRepeatLength) && (actualRepeatLength <= maxRepeatLength) )
            {   if (CRISPRUtil.hasNonRepeatingSpacers(candidateCRISPR, spacerToSpacerMaxSimilarity))
               {
                  if (CRISPRUtil.hasSimilarlySizedSpacers(candidateCRISPR, spacerToSpacerLengthDiff, spacerToRepeatLengthDiff))
                  {
                     candidateCRISPR = CRISPRUtil.checkFlank("left", candidateCRISPR, minSpacerLength, 30, spacerToSpacerMaxSimilarity, .70);
                     candidateCRISPR = CRISPRUtil.checkFlank("right", candidateCRISPR, minSpacerLength, 30, spacerToSpacerMaxSimilarity, .70);

                     candidateCRISPR = CRISPRUtil.trim(candidateCRISPR, minRepeatLength);

                     CRISPRVector.addElement(candidateCRISPR);
                     repeatsFound = true;

                     //we may skip current CRISPR (assuming CRISPRs are not interleaved)
                     j = candidateCRISPR.end() + 1;
                  }

               }
            }
         }
      }

      long repeatSearchEnd = System.currentTimeMillis();
      //System.out.println("Time to search for repeats:  " + (repeatSearchEnd - repeatSearchStart) + " ms");
      //System.out.println(CRISPRVector.size() + " possible CRISPR(s) found" + "\n");


      // ********************** Display CRISPR elements ********************** //
      // ********************************************************************* //
      try
      {

         FileOutputStream outputFileStream;
         PrintStream out;

         FileOutputStream outputGffFileStream;
         PrintStream gffOut = null;

         if (screenDisplay == 1)
            out = System.out;
         else
         {
            if ( outputFileName == "" )
               outputFileName = "a.out";


            File outputFile = new File(outputFileName);
            if ( readNum == 1 && outputFile.exists() )
            {
               boolean success = outputFile.delete();
               if (!success)
                  throw new IllegalArgumentException("Error: Could not delete file '" + outputFile + "'");
            }

            outputFileStream = new FileOutputStream(outputFile, true);
            out = new PrintStream(outputFileStream);


            if (! outputGffFileName.equals(""))
            {
              File outputGffFile = new File(outputGffFileName);
              if ( readNum == 1 && outputGffFile.exists() )
              {
                 boolean success = outputFile.delete();
                 if (!success)
                    throw new IllegalArgumentException("Error: Could not delete file '" + outputFile + "'");
              }

              outputGffFileStream = new FileOutputStream(outputGffFile, true);
              gffOut = new PrintStream(outputGffFileStream);
              gffOut.println("##gff-version 3");
            	printGffHeader = false;
            }
         }

         if (repeatsFound)
         {
             if(outputformat == 0) {
                 out.print("Sequence '" + sequence.getName() + "' (" + sequence.length() + " bp)\n");
                 out.print("\n");
             } else {
            	 if(printGffHeader) {
            		 out.println("##gff-version 3");
            		 printGffHeader = false;
            	 }
             }
        	 //int repeatLength, numRepeats, numSpacers;
            CRISPR currCRISPR;

            //String repeat, spacer, prevSpacer;
            //repeat = spacer = prevSpacer = "";

            //add 1 to each position, to offset programming languagues that begin at 0 rather than 1
            for (int k = 0; k < CRISPRVector.size(); k++)
            {
                currCRISPR = (CRISPR)CRISPRVector.elementAt(k);
                totalCrisprCount++;
                if(outputformat > 0 && gffOut == null ) {
                  printGff(out, sequence, currCRISPR);
                } else {
                  printTable(out, currCRISPR);
                }

                if (gffOut != null) {
                  printGff(gffOut, sequence, currCRISPR);
                }
                if(printSpacers) {
                	for (int i = 0; i < currCRISPR.numSpacers(); ++i) {
                    	spacers.print(">"+ sequence.getName() +"_CRISPR_"+ (totalCrisprCount) +"_spacer_"+ (i+1)+"\n" + currCRISPR.spacerStringAt(i) + "\n");
                	}
                }


            }
            if(outputformat == 0) {
            	out.print("Time to find repeats: " + (repeatSearchEnd - repeatSearchStart) + " ms\n\n");
            	out.print("\n");
            }
         }




         if (screenDisplay == 0)
            out.close();

      }
      catch (Exception e)   {   
        System.err.println ("--Error writing to file-- \n");   
        e.printStackTrace(System.err);
      }

      return true;
   }

   private boolean printGff(PrintStream out, DNASequence sequence, CRISPR currCRISPR) {
     String crispr_id = "CRISPR" + totalCrisprCount;
     out.print(sequence.getName() + "\tminced:" + minced.VERSION + "\trepeat_region\t");
     out.print((currCRISPR.start() + 1) + "\t" + (currCRISPR.end() + 1) + "\t");
     out.print(currCRISPR.numRepeats() + "\t.\t.\tID="+ crispr_id + ";rpt_type=direct;rpt_family=CRISPR;rpt_unit_seq="+ currCRISPR.repeatStringAt(1));
     out.print("\n");
     if(outputformat == 2) {
       out.print(currCRISPR.toGff(sequence.getName(), crispr_id));
     }
    return true;
   }

   private boolean printTable(PrintStream out, CRISPR currCRISPR) {
     out.print("CRISPR " + totalCrisprCount + "   Range: " + (currCRISPR.start() + 1) + " - " +  (currCRISPR.end() + 1) + "\n");
     out.print(currCRISPR.toString());
     out.print("Repeats: " + currCRISPR.numRepeats() + "\t" +  "Average Length: " + currCRISPR.averageRepeatLength() + "\t\t");
     out.print("Average Length: " +  currCRISPR.averageSpacerLength() + "\n\n");
     return true;
   }

}

