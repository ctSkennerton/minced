import java.io.*;
import java.io.File;
import java.util.Vector;

public class CRISPRFinder
{
   private String inputFileName;
   private String outputFileName;

   private int screenDisplay;
   private int minNumRepeats;
   private int minRepeatLength;
   private int maxRepeatLength;
   private int minSpacerLength;
   private int maxSpacerLength;
   private int searchWindowLength;

   DNASequence sequence = null;
   int sequenceLength = 0;


   public CRISPRFinder(String _inputFileName, String _outputFileName, int _screenDisplay, int _minNumRepeats, int _minRepeatLength, int _maxRepeatLength, int _minSpacerLength, int _maxSpacerLength, int _searchWindowLength)
   {
      inputFileName = _inputFileName;
      outputFileName = _outputFileName;

      screenDisplay = _screenDisplay;
      minNumRepeats = _minNumRepeats;
      minRepeatLength = _minRepeatLength;
      maxRepeatLength = _maxRepeatLength;
      minSpacerLength = _minSpacerLength;
      maxSpacerLength = _maxSpacerLength;
      searchWindowLength = _searchWindowLength;
   }


   public boolean goCRISPRFinder()
   {
      File inputFile = new File(inputFileName);
      if (! inputFile.exists())
      {
         System.out.println("File name does not exist:  " + inputFile.getPath());
         return false;
      }
      System.out.println("\n\nReading file '" + inputFile.getPath() + "'");
      FASTAReader fastaReader = new FASTAReader(inputFile.getPath());
      if (! fastaReader.isFASTA())
      {
         System.out.println("File '" + inputFile.getPath() + "' does not seem to be FASTA-formatted.");
         return false;
      }

      while ( fastaReader.read() )
      {
         sequence = new DNASequence(fastaReader.getSequence(), fastaReader.getName(), fastaReader.getDesc());
         System.out.println("Sequence '" + sequence.getName() + "' (" + sequence.length() + " bp)");

         /*
         // Check that sequence is a valid DNA sequence
         if ( ! sequence.isDNASequence() )
         {
            System.out.println(sequence.getErrorLog() + "Skipping this sequence!");
            System.out.println("");
            continue;
         }
         */

         try
         {
            findRepeats(sequence, fastaReader.getNum());
         }
         catch (Exception e)
         {
            System.out.println ("Error processing input file '" + inputFile.getPath() + "'. Please, check contents.\n");
         }

      }

      System.out.println("Done!");
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
         System.out.println("Changing window length to " + searchWindowLength + " instead of " + oldSearchWindowlength);
      }

      double spacerToSpacerMaxSimilarity = 0.62;
      int spacerToSpacerLengthDiff = 12;
      int spacerToRepeatLengthDiff = 30;

      //the mumber of bases that can be skipped while we still guarantee that the entire search
      //window will at some point in its iteration thru the sequence will not miss a any repeat
      int skips = minRepeatLength - (2 * searchWindowLength - 1);
      if (skips < 1)
         skips = 1;

      System.out.println("Searching for repeats...");
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
            candidateCRISPR = CRISPRUtil.getActualRepeatLength(candidateCRISPR, searchWindowLength, minSpacerLength);
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
      System.out.println("Time to search for repeats:  " + (repeatSearchEnd - repeatSearchStart) + " ms");
      System.out.println(CRISPRVector.size() + " possible CRISPR(s) found" + "\n");


      // ********************** Display CRISPR elements ********************** //
      // ********************************************************************* //
      try
      {

         FileOutputStream outputFileStream;
         PrintStream out;
         
         if (screenDisplay == 1)
            out = System.out;
         else
         {
            if ( outputFileName.equals("") )
               outputFileName = "a.out";

            System.out.println("Writing results in file '" + outputFileName + "'");
            System.out.println("");

            File outputFile = new File(outputFileName);
            if ( readNum == 1 && outputFile.exists() )
            {
               boolean success = outputFile.delete();
               if (!success)
                  throw new IllegalArgumentException("Error: Could not delete file '" + outputFile + "'");
            }

            outputFileStream = new FileOutputStream(outputFile, true);
            out = new PrintStream(outputFileStream);
         }

         out.print("Sequence '" + sequence.getName() + "' (" + sequence.length() + " bp)\n");
         out.print("\n");

         if (repeatsFound)
         {
            int repeatLength, numRepeats, numSpacers;
            CRISPR currCRISPR;

            String repeat, spacer, prevSpacer;
            repeat = spacer = prevSpacer = "";

            //add 1 to each position, to offset programming languagues that begin at 0 rather than 1
            for (int k = 0; k < CRISPRVector.size(); k++)
            {
               currCRISPR = (CRISPR)CRISPRVector.elementAt(k);
               out.print("CRISPR " + (k + 1) + "   Range: " + (currCRISPR.start() + 1) + " - " +  (currCRISPR.end() + 1) + "\n");
               out.print(currCRISPR.toString());
               out.print("Repeats: " + currCRISPR.numRepeats() + "\t" +  "Average Length: " + currCRISPR.averageRepeatLength() + "\t\t");
               out.print("Average Length: " +  currCRISPR.averageSpacerLength() + "\n\n");
               out.print("\n\n");
            }
            //out.print("Time to find repeats: " + (repeatSearchEnd - repeatSearchStart) + " ms\n\n");
         }

         out.print("\n");

         if (!repeatsFound)
            out.print("No CRISPR elements were found.\n\n");

         if (screenDisplay == 0)
            out.close();

      }
      catch (Exception e)   {   System.err.println ("--Error writing to file-- \n");   }

      return true;
   }
}
