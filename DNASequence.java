import java.io.*;
import java.util.*;
import java.util.regex.*;

public class DNASequence
{
   public static String seq;
   private String name;
   private String desc;
   private String errorLog = "";
   private Pattern notDnaRe = Pattern.compile("([^ACGT])");
   private IntervalSearchTree mask = new IntervalSearchTree();

   public DNASequence(String sequence)
   {
      seq  = sequence.toUpperCase();
      name = "";
      desc = "";
   }

   public DNASequence(String sequence, String _name)
   {
      seq  = sequence.toUpperCase();
      name = _name;
      desc = "";
   }

   public DNASequence(String sequence, String _name, String _desc)
   {
      seq  = sequence.toUpperCase();
      name = _name;
      desc = _desc;
   }

   public int length()
   {
      return seq.length();
   }

   public String getName()
   {
      return name;
   }

   public String getDesc()
   {
      return desc;
   }

   public String toString()
   {
      return seq;
   }

   public boolean isDNASequence()
   {
      int length = seq.length();

      if (length == 0)
      {
         errorLog = "Not a valid DNA sequence. Sequence contains no characters.";
         return false;
      }

      Matcher m = notDnaRe.matcher(seq);
      boolean isNotDNA = m.find();
      if ( isNotDNA )
      {
         errorLog = "Not a valid DNA sequence. Invalid character '" + m.group() + "' found at position " + m.start() + ". ";
         return false;
      }
      else
      {
         return true;
      }
   }

   /*
   //method has not been tested
   //may contain alphabetic characters other than acgt
   public boolean isSequence()
   {
      int i = 0;
      int length = seq.length();

      if (length == 0)
      {   errorLog = "Not a valid sequence.  Sequence contains no characters.";
         return false;
      }

      while (i < length)
      {
         if (Character.isLetter(seq.charAt(i)))
         {   System.out.println(seq.charAt(i) + " " + Character.isLetter(seq.charAt(i)));
            i++;
         }
         else
         {   errorLog = "Not a valid sequence.  Error found at position " + i + " scanning a non-alphabetic character " + seq.charAt(i);
            return false;
         }
      }
      return true;
   }
   */

   public String getErrorLog()
   {   return errorLog;
   }


   public String substring(int beginIndex, int endIndex)
   {
      return seq.substring(beginIndex, endIndex);
   }

   public void toUpperCase()
   {
      seq = seq.toUpperCase();
   }

   public static int getHammingDistance(String seq1, String seq2)
   {
      int length = seq1.length();
      int hammingDistance = 0;

      if (seq1.length() != seq2.length())
      {
         length = min(seq1.length(), seq2.length());
         hammingDistance = Math.abs(seq1.length() - seq2.length());
      }

      for (int i =0; i < length; i++)
      {
         if ( seq1.charAt(i) != seq2.charAt(i))
         {
            hammingDistance++;
         }
      }

      return hammingDistance;
   }


   //*****************************
   // Compute Levenshtein distance  Michael Gilleland, Merriam Park Software www.merriampark.com/ld.htm
   //*****************************
   public static int getLevenshteinDistance1 (String s, String t)
   {
      int d[][]; // matrix
      int n;     // length of s
      int m;     // length of t
      int i;     // iterates through s
      int j;     // iterates through t
      char s_i;  // ith character of s
      char t_j;  // jth character of t
      int cost;  // cost

      // Step 1
      n = s.length ();
      m = t.length ();
      if (n == 0) {
         return m;
      }
      if (m == 0) {
         return n;
      }
      d = new int[n+1][m+1];

      // Step 2
      for (i = 0; i <= n; i++) {
         d[i][0] = i;
      }

      for (j = 0; j <= m; j++) {
         d[0][j] = j;
      }

      // Step 3
      for (i = 1; i <= n; i++) {

         s_i = s.charAt (i - 1);

         // Step 4
         for (j = 1; j <= m; j++) {

            t_j = t.charAt (j - 1);

            // Step 5
            if (s_i == t_j) {
               cost = 0;
            }
            else {
               cost = 1;
            }

            // Step 6
            d[i][j] = Minimum (d[i-1][j]+1, d[i][j-1]+1, d[i-1][j-1] + cost);

         }

      }

      // Step 7
      return d[n][m];

   }


   //Chas Emerick  http://www.merriampark.com/ldjava.htm
   public static int getLevenshteinDistance2 (String s, String t) {
      if (s == null || t == null)
      {
         throw new IllegalArgumentException("Strings must not be null");
      }

      /*
      The difference between this impl. and the previous is that, rather
      than creating and retaining a matrix of size s.length()+1 by t.length()+1,
      we maintain two single-dimensional arrays of length s.length()+1.  The first, d,
      is the 'current working' distance array that maintains the newest distance cost
      counts as we iterate through the characters of String s.  Each time we increment
      the index of String t we are comparing, d is copied to p, the second int[].  Doing so
      allows us to retain the previous cost counts as required by the algorithm (taking
      the minimum of the cost count to the left, up one, and diagonally up and to the left
      of the current cost count being calculated).  (Note that the arrays aren't really
      copied anymore, just switched...this is clearly much better than cloning an array
      or doing a System.arraycopy() each time through the outer loop.)

      Effectively, the difference between the two implementations is this one does not
      cause an out of memory condition when calculating the LD over two very large strings.
      */

      int n = s.length(); // length of s
      int m = t.length(); // length of t

      if (n == 0) {
         return m;
      } else if (m == 0) {
         return n;
      }

      int p[] = new int[n+1]; //'previous' cost array, horizontally
      int d[] = new int[n+1]; // cost array, horizontally
      int _d[]; //placeholder to assist in swapping p and d

      // indexes into strings s and t
      int i; // iterates through s
      int j; // iterates through t

      char t_j; // jth character of t

      int cost; // cost

      for (i = 0; i<=n; i++) {
         p[i] = i;
      }

      for (j = 1; j<=m; j++) {
         t_j = t.charAt(j-1);
         d[0] = j;

         for (i=1; i<=n; i++) {
            cost = s.charAt(i-1)==t_j ? 0 : 1;
            // minimum of cell to the left+1, to the top+1, diagonally left and up +cost
            d[i] = Math.min(Math.min(d[i-1]+1, p[i]+1),  p[i-1]+cost);
         }

         // copy current distance counts to 'previous row' distance counts
         _d = p;
         p = d;
         d = _d;
      }

      // our last action in the above loop was to switch d and p, so p now
      // actually has the most recent cost counts
      return p[n];
   }

   public static double getSimilarity(String s1, String s2)
   {   int maxLength = max(s1.length(), s2.length());
      double similarity = 1.0 - (double)getLevenshteinDistance1(s1, s2)/maxLength;
      return similarity;
   }

   private static boolean patternMatches(String pattern1, String pattern2, double confidence)
   {   double patternSimilarity = getSimilarity(pattern1, pattern2);
      if (patternSimilarity >= confidence)
         return true;
      else
         return false;
   }


   private static int min (int[] array)
   {   int min = array[0];
      int minIndex = 0;

      for (int i = 0; i < array.length; i++)
      {   if (array[i] < min)
         {   min = array[i];
            minIndex = i;
         }
      }
      return minIndex;
   }

   private static int max (int[] array)
   {   int max = array[0];
      int maxIndex = 0;

      for (int i = 0; i < array.length; i++)
      {   if (array[i] > max)
         {   max = array[i];
            maxIndex = i;
         }
      }
      return maxIndex;
   }

   private static int min(int n1, int n2)
   {   if (n1 < n2)
         return n1;
      else
         return n2;
   }

   private static int max(int n1, int n2)
   {   if (n1 < n2)
         return n2;
      else
         return n1;
   }

   private static int Minimum (int a, int b, int c)
   {
      int mi;

      mi = a;
      if (b < mi)
         mi = b;

      if (c < mi)
         mi = c;

      return mi;

   }

   public void mask(int minimum)
   {
       for(int i = 0; i < seq.length(); i++)
       {
           int n = 0;
           int j = i + 1;
           while(j < seq.length() && seq.charAt(i) == seq.charAt(j))
           {
               ++n;
               ++j;
           }
           if(n >= minimum)
           {
               // add in the interval and save the end location as the data
               mask.add(i, j);
           }
           i = j;

       }
   }

   public int masked(int start, int end)
   {
       return mask.overlapEnd(start, end);
   }

}
