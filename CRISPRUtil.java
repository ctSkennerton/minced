import java.util.*;

public class CRISPRUtil
{
   //identified repeats may represent only a subset of a larger repeat.  this method extends these
   //repeats as long as they continue to match within some range.  assumes there are at least two repeats
   public static CRISPR getActualRepeatLength(CRISPR candidateCRISPR, int searchWindowLength, int minSpacerLength, int minNumRepeats)
   {
      int firstRepeatStartIndex = candidateCRISPR.repeatAt(0);
      int lastRepeatStartIndex = candidateCRISPR.repeatAt(candidateCRISPR.numRepeats()-1);

      int shortestRepeatSpacing = candidateCRISPR.repeatAt(1) - candidateCRISPR.repeatAt(0);
      for (int i = 0; i < candidateCRISPR.numRepeats() - 1; i++)
      {
         int currRepeatIndex = candidateCRISPR.repeatAt(i);
         int nextRepeatIndex = candidateCRISPR.repeatAt(i + 1);
         int currRepeatSpacing = nextRepeatIndex - currRepeatIndex;
         if (currRepeatSpacing < shortestRepeatSpacing)
            shortestRepeatSpacing = currRepeatSpacing;
      }

      int sequenceLength = DNASequence.seq.length();

      int rightExtensionLength = searchWindowLength;
      int maxRightExtensionLength = shortestRepeatSpacing - minSpacerLength;


      int currRepeatStartIndex;
      String currRepeat;
      int charCountA, charCountC, charCountT, charCountG;
      charCountA = charCountC = charCountT = charCountG = 0;
      double threshold;
      boolean done = false;

      threshold = .75;

      //(from the right side) extend the length of the repeat to the right as long as the last base of all repeats are at least threshold
      while (!done && (rightExtensionLength <= maxRightExtensionLength))
      {
         if (lastRepeatStartIndex + rightExtensionLength >= sequenceLength)
         {
            if (candidateCRISPR.numRepeats() - 1 > minNumRepeats)
            {
               candidateCRISPR.removeRepeat(candidateCRISPR.lastRepeat());
               lastRepeatStartIndex = candidateCRISPR.repeatAt(candidateCRISPR.numRepeats()-1);
            }
            else
            {
               done = true;
               break;
            }
         }
         for (int k = 0; k < candidateCRISPR.numRepeats(); k++ )
         {
            currRepeatStartIndex = candidateCRISPR.repeatAt(k);
            //System.out.println("currRepeatStartIndex: " + currRepeatStartIndex + " currRepeatStartIndex + rightExtensionLength " + (currRepeatStartIndex + rightExtensionLength) + " maxRightExtensionLength "+maxRightExtensionLength+" sequence length "+sequenceLength);
            currRepeat = DNASequence.seq.substring(currRepeatStartIndex, currRepeatStartIndex + rightExtensionLength);
            char lastChar = currRepeat.charAt(currRepeat.length() - 1);

            if (lastChar == 'A')   charCountA++;
            if (lastChar == 'C')   charCountC++;
            if (lastChar == 'T')   charCountT++;
            if (lastChar == 'G')   charCountG++;
         }

         double percentA = (double)charCountA/candidateCRISPR.numRepeats();
         double percentC = (double)charCountC/candidateCRISPR.numRepeats();
         double percentT = (double)charCountT/candidateCRISPR.numRepeats();
         double percentG = (double)charCountG/candidateCRISPR.numRepeats();

         if ( (percentA >= threshold) || (percentC >= threshold) || (percentT >= threshold) || (percentG >= threshold) )
         {
            rightExtensionLength++;
            charCountA = charCountC = charCountT = charCountG = 0;
         }
         else
         {
            done = true;
         }
      }
      rightExtensionLength--;


      int leftExtensionLength = 0;
      charCountA = charCountC = charCountT = charCountG = 0;
      done = false;

      int maxLeftExtensionLength = shortestRepeatSpacing - minSpacerLength - rightExtensionLength;

      //(from the left side) extends the length of the repeat to the left as long as the first base of all repeats is at least threshold
      while (!done && (leftExtensionLength <= maxLeftExtensionLength) )
      {
         if (firstRepeatStartIndex - leftExtensionLength < 0)
         {
            if (candidateCRISPR.numRepeats() - 1 > minNumRepeats)
            {
               candidateCRISPR.removeRepeat(candidateCRISPR.firstRepeat());
               firstRepeatStartIndex = candidateCRISPR.repeatAt(0);
            }
            else
            {
               done = true;
               break;
            }
         }
         for (int k = 0; k < candidateCRISPR.numRepeats(); k++ )
         {
            currRepeatStartIndex = candidateCRISPR.repeatAt(k);
            char firstChar = DNASequence.seq.charAt(currRepeatStartIndex - leftExtensionLength);

            if (firstChar == 'A')    charCountA++;
            if (firstChar == 'C')    charCountC++;
            if (firstChar == 'T')    charCountT++;
            if (firstChar == 'G')    charCountG++;
         }

         double percentA = (double)charCountA/candidateCRISPR.numRepeats();
         double percentC = (double)charCountC/candidateCRISPR.numRepeats();
         double percentT = (double)charCountT/candidateCRISPR.numRepeats();
         double percentG = (double)charCountG/candidateCRISPR.numRepeats();

         if ( (percentA >= threshold) || (percentC >= threshold) || (percentT >= threshold) || (percentG >= threshold) )
         {
            leftExtensionLength++;
            charCountA = charCountC = charCountT = charCountG = 0;
         }
         else
         {
            done = true;
         }
      }
      leftExtensionLength--;

      // Ok to suppress warnings since we know that we are casting integers
      @SuppressWarnings("unchecked")
      Vector<Integer> newPositions = (Vector<Integer>)(candidateCRISPR.repeats()).clone();

      for (int m = 0; m < newPositions.size(); m++)
      {
         int newValue = candidateCRISPR.repeatAt(m) - leftExtensionLength;
         newPositions.setElementAt(new Integer(newValue), m);
      }

      int actualPatternLength = rightExtensionLength + leftExtensionLength;


      return new CRISPR(newPositions, actualPatternLength);

   }


   public static CRISPR trim(CRISPR candidateCRISPR, int minRepeatLength)
   {
      int numRepeats = candidateCRISPR.numRepeats();
      int left = candidateCRISPR.start();
      int right = candidateCRISPR.end();

      String currRepeat;
      int charCountA, charCountC, charCountT, charCountG;
      charCountA = charCountC = charCountT = charCountG = 0;
      boolean done = false;

      //trim from right
      while (!done && (candidateCRISPR.repeatLength() > minRepeatLength) )
      {
         for (int k = 0; k < candidateCRISPR.numRepeats(); k++ )
         {
            currRepeat = candidateCRISPR.repeatStringAt(k);
            char lastChar = currRepeat.charAt(currRepeat.length() - 1);

            if (lastChar == 'A')   charCountA++;
            if (lastChar == 'C')   charCountC++;
            if (lastChar == 'T')   charCountT++;
            if (lastChar == 'G')   charCountG++;
         }

         double percentA = (double)charCountA/candidateCRISPR.numRepeats();
         double percentC = (double)charCountC/candidateCRISPR.numRepeats();
         double percentT = (double)charCountT/candidateCRISPR.numRepeats();
         double percentG = (double)charCountG/candidateCRISPR.numRepeats();

         if ( (percentA < .75) && (percentC < .75) && (percentT < .75) && (percentG < .75) )
         {
            candidateCRISPR.setRepeatLength(candidateCRISPR.repeatLength() - 1);
            charCountA = charCountC = charCountT = charCountG = 0;
         }
         else
         {
            done = true;
         }
      }



      charCountA = charCountC = charCountT = charCountG = 0;
      done = false;

      //trim from left
      while (!done && (candidateCRISPR.repeatLength() > minRepeatLength) )
      {
         for (int k = 0; k < candidateCRISPR.numRepeats(); k++ )
         {
            currRepeat = candidateCRISPR.repeatStringAt(k);
            char firstChar = currRepeat.charAt(0);

            if (firstChar == 'A')   charCountA++;
            if (firstChar == 'C')   charCountC++;
            if (firstChar == 'T')   charCountT++;
            if (firstChar == 'G')   charCountG++;
         }

         double percentA = (double)charCountA/candidateCRISPR.numRepeats();
         double percentC = (double)charCountC/candidateCRISPR.numRepeats();
         double percentT = (double)charCountT/candidateCRISPR.numRepeats();
         double percentG = (double)charCountG/candidateCRISPR.numRepeats();

         if ( (percentA < .75) && (percentC < .75) && (percentT < .75) && (percentG < .75) )
         {
            for (int m = 0; m < candidateCRISPR.numRepeats(); m++ )
            {
               int newValue = candidateCRISPR.repeatAt(m) + 1;
               candidateCRISPR.setRepeatAt(newValue, m);
            }
            candidateCRISPR.setRepeatLength(candidateCRISPR.repeatLength() - 1);
            charCountA = charCountC = charCountT = charCountG = 0;
         }
         else
         {
            done = true;
         }
      }

      return candidateCRISPR;
   }


   public static boolean hasSimilarlySizedSpacers(CRISPR candidateCRISPR, int spacerToSpacerLengthOffset, int spacerToRepeatLengthOffset)
   {
      int initialSpacerLength = candidateCRISPR.spacerStringAt(0).length();
      int repeatLength = candidateCRISPR.repeatLength();

      for (int i = 0 ; i < candidateCRISPR.numSpacers(); i++)
      {
         int currSpacerLength = candidateCRISPR.spacerStringAt(i).length();
         //checks that each spacer is of similar size to other spacers
         if ( Math.abs(currSpacerLength - initialSpacerLength) >  spacerToSpacerLengthOffset )
         {
            return false;
         }

         //checks that each spacer is of similar size to the repeats
         if ( Math.abs(currSpacerLength - repeatLength) > spacerToRepeatLengthOffset)
         {
            return false;
         }

      }
      return true;
   }


   //checks first five spacers
   public static boolean hasNonRepeatingSpacers(CRISPR candidateCRISPR, double maxSimilarity)
   {
      //assumes at least two elements
      String firstRepeat = candidateCRISPR.repeatStringAt(0);
      String firstSpacer = candidateCRISPR.spacerStringAt(0);

      if (candidateCRISPR.numRepeats() >= 3)
      {
         int i = 0;
         while ( (i < candidateCRISPR.numSpacers() - 1) )
         {
            if (i == 4)  //only check first 5 spacers
               return true;

            String currSpacer = candidateCRISPR.spacerStringAt(i);
            String nextSpacer = candidateCRISPR.spacerStringAt(i + 1);
            String currRepeat = candidateCRISPR.repeatStringAt(i);

            //spacers should be different
            if ( DNASequence.getSimilarity(currSpacer, nextSpacer) > maxSimilarity )
               return false;

            //repeats should also be different from spacers, otherwise may be tandem repeat
            if ( DNASequence.getSimilarity(currRepeat, currSpacer) > maxSimilarity )
               return false;

            i++;
         }

         //checks last repeat/spacer
         if ( DNASequence.getSimilarity(candidateCRISPR.repeatStringAt(i), candidateCRISPR.spacerStringAt(i)) > maxSimilarity )
            return false;

         return true;
      }

      //we check that the spacer is different from the repeat
      else if (candidateCRISPR.numRepeats() == 2)
      {
         if (firstSpacer.equals(""))
            return false;
         else
            return ( DNASequence.getSimilarity(firstSpacer, firstRepeat) < maxSimilarity );
      }

      else
      {
         return false;
      }

   }


   public static CRISPR checkFlank(String side, CRISPR candidateCRISPR, int minSpacerLength, int scanRange, double spacerToSpacerMaxSimilarity, double confidence)
   {
      boolean moreToSearch = true;

      while (moreToSearch)
      {
         int result = scan(side, candidateCRISPR, minSpacerLength, scanRange, confidence);
         if (result > 0)  //if another repeat found on flank
         {
            if (side.equals("left"))
               candidateCRISPR.insertRepeatAt(result, 0);
            else if (side.equals("right"))
               candidateCRISPR.addRepeat(result);
         }
         else
            moreToSearch = false;
      }

      return candidateCRISPR;
   }

   /*
      scan to the right and left of the first and last repeat to see if there is a region
      that is similar to the repeats.  necessary in case we missed a repeat because of
      inexact matches or a result of one of the filters
   */
   private static int scan(String side, CRISPR candidateCRISPR, int minSpacerLength, int scanRange, double confidence)
   {
      int repeatSpacing1, repeatSpacing2, avgRepeatSpacing;
      int firstRepeatIndex, lastRepeatIndex, candidateRepeatIndex;
      String repeatString, candidateRepeatString, newCandidateRepeatString;

      int repeatLength = candidateCRISPR.repeatLength();
      int numRepeats = candidateCRISPR.numRepeats();
      int sequenceLength = DNASequence.seq.length();

      firstRepeatIndex = candidateCRISPR.repeatAt(0);
      lastRepeatIndex = candidateCRISPR.repeatAt(numRepeats-1);

      if (side.equals("left"))
      {
         repeatString = candidateCRISPR.repeatStringAt(0);
         repeatSpacing1 = candidateCRISPR.repeatSpacing(0, 1);
         if (numRepeats >= 3)
         {
            repeatSpacing2 = candidateCRISPR.repeatSpacing(1, 2);
            avgRepeatSpacing = (repeatSpacing1 + repeatSpacing2)/2;
         }
         else
            avgRepeatSpacing = repeatSpacing1;

         candidateRepeatIndex = firstRepeatIndex - avgRepeatSpacing;
      }

      else //if (side.equals("right"))
      {
         repeatString = candidateCRISPR.repeatStringAt(numRepeats-1);
         repeatSpacing1 = candidateCRISPR.repeatSpacing(numRepeats-2, numRepeats-1);
         if (numRepeats >= 3)
         {
            repeatSpacing2 = candidateCRISPR.repeatSpacing(numRepeats-3, numRepeats-2);
            avgRepeatSpacing = (repeatSpacing1 + repeatSpacing2)/2;
         }
         else
            avgRepeatSpacing = repeatSpacing1;

         candidateRepeatIndex = lastRepeatIndex + avgRepeatSpacing;
      }

      int begin = candidateRepeatIndex - scanRange;
      int end   = candidateRepeatIndex + scanRange;

      /******************** range checks ********************/
      //check that we do not search too far within an existing repeat when scanning right and left
      int scanLeftMaxEnd    = firstRepeatIndex - repeatLength - minSpacerLength;
      int scanRightMinBegin = lastRepeatIndex + repeatLength + minSpacerLength;

      if (side.equals("left"))
      {
         if (end > scanLeftMaxEnd)
            end = scanLeftMaxEnd;
      }

      if (side.equals("right"))
      {
         if (begin < scanRightMinBegin)
            begin = scanRightMinBegin;
      }

      //out of bounds check for scanning left
      if ( (begin) < 0)
         return 0;

      //out of bounds check for scanning right
      if ( (begin + repeatLength) > sequenceLength)
         return 0;
      if ( (end + repeatLength) > sequenceLength)
         end = sequenceLength - repeatLength;

      if ( begin >= end)
         return 0;
      /******************** end range checks ********************/

      int[] array = new int[end - begin + 1];

      int index = 0;
      for (int i = begin; i <= end; i++)
      {
         candidateRepeatString = DNASequence.seq.substring(i, i + repeatLength);
         array[index] = DNASequence.getHammingDistance(repeatString, candidateRepeatString);
         index++;
      }

      //min(array) returns the index of the smallest value in array  in this case, it refers to
      //the candidate string theat is closest to the repeatString.  uses hamming distance as levenshteinDistance is not useful for this particular task
      int newCandidateRepeatIndex = begin + min(array);
      newCandidateRepeatString = DNASequence.seq.substring(newCandidateRepeatIndex, newCandidateRepeatIndex + repeatLength);

      boolean match = patternMatches(repeatString, newCandidateRepeatString, confidence);

      if (match)
         return newCandidateRepeatIndex;
      else
         return 0;

   }

   public static CRISPR scanRight(CRISPR candidateCRISPR, String pattern, int minSpacerLength, int scanRange, SearchUtil searchUtil)
   {
      int numRepeats = candidateCRISPR.numRepeats();
      int patternLength = pattern.length();
      int sequenceLength = DNASequence.seq.length();

      int lastRepeatIndex = candidateCRISPR.repeatAt(numRepeats-1);
      int secondToLastRepeatIndex = candidateCRISPR.repeatAt(numRepeats-2);
      int repeatSpacing = lastRepeatIndex - secondToLastRepeatIndex;

      int candidateRepeatIndex, beginSearch, endSearch, position;
      String text = "";

      boolean moreToSearch = true;
      while (moreToSearch)
      {
         candidateRepeatIndex = lastRepeatIndex + repeatSpacing;
         beginSearch = candidateRepeatIndex - scanRange;
         endSearch = candidateRepeatIndex + patternLength + scanRange;

         /******************** range checks ********************/
         //check that we do not search too far within an existing repeat when scanning right
         int scanRightMinBegin = lastRepeatIndex + patternLength + minSpacerLength;

         if (beginSearch < scanRightMinBegin)
            beginSearch = scanRightMinBegin;

         //System.out.print("beginSearch " + beginSearch + "  " + "endSearch" + endSearch);
         if (beginSearch > sequenceLength - 1)
            return candidateCRISPR;
         if (endSearch > sequenceLength)
            endSearch = sequenceLength;

         if ( beginSearch >= endSearch)
            return candidateCRISPR;
         /******************** end range checks ********************/

         text = DNASequence.seq.substring(beginSearch, endSearch);
         position = searchUtil.boyer_mooreSearch(text, pattern);

         if (position >= 0)
         {
            candidateCRISPR.addRepeat(beginSearch + position);
            secondToLastRepeatIndex = lastRepeatIndex;
            lastRepeatIndex = beginSearch + position;
            repeatSpacing = lastRepeatIndex - secondToLastRepeatIndex;
            if (repeatSpacing < (minSpacerLength + patternLength))
               moreToSearch = false;
         }
         else
            moreToSearch = false;
      }

      return candidateCRISPR;
   }


   private static boolean patternMatches(String pattern1, String pattern2, double confidence)
   {
      double patternSimilarity = DNASequence.getSimilarity(pattern1, pattern2);
      if (patternSimilarity >= confidence)
         return true;
      else
         return false;
   }

   private static int min (int[] array)
   {
      int min = array[0];
      int minIndex = 0;

      for (int i = 0; i < array.length; i++)
      {
         if (array[i] < min)
         {
            min = array[i];
            minIndex = i;
         }
      }
      return minIndex;
   }
}

