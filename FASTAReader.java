import java.io.*;
import java.util.*;

public class FASTAReader
{
   private String fileName;
   private String sequence;
   private String name;
   private String desc;
   private String newHeader;
   private int readNum;
   private BufferedReader inputFile;

   public FASTAReader(String f)
   {
      fileName  = f;
      readNum   = 0;
      sequence  = "";
      name      = "";
      desc      = "";
      newHeader = "";
   }

   //does not check that all characters are ACGT
   public boolean isFASTA()
   {
      try
      {
         BufferedReader testFile = new BufferedReader(new FileReader(fileName));  

         String firstLine = testFile.readLine();
         while (firstLine.length() == 0)
         {  
            firstLine = testFile.readLine();
         }
         
         if (firstLine == null)
         {
            testFile.close();
            return false;
         }
         if (firstLine.charAt(0) != '>')
         {
            testFile.close();
            return false;
         }

         String secondLine = testFile.readLine();
         if (secondLine == null)
         {
            testFile.close();
            return false;
         }

         testFile.close();
      } catch (Exception e) { e.printStackTrace(); }

      return true;
   }

   public boolean read()
   {
      try
      {
         if ( readNum == 0 )
            inputFile = new BufferedReader(new FileReader(fileName));
         
         readNum++;
         sequence  = "";

         StringBuffer sb = new StringBuffer();

         while( true )
         {
            String currLine = inputFile.readLine();

            // Process informative line
            if ( currLine == null )
            {
               if ( sb.toString().length() > 0 )
               {
                  // Last sequence
                  sequence = sb.toString();
                  return true;
               }
               else
               {
                  // No more sequences
                  inputFile.close();
                  return false;
               }
            }
            else
            {
               // Remove leading and trailing spaces
               currLine = currLine.trim();

               if ( currLine.length() == 0 )
               {
                  // Skip empty lines
                  continue;
               }
               else if ( currLine.charAt(0) == '>' )
               {               
                  // Prepare new sequence
                  newHeader = currLine;
                  // Finish previous sequence
                  if ( sb.toString().length() > 0 ) 
                  {
                     sequence = sb.toString();
                     return true;
                  }
               }
               else
               {
                  // Continuation of a sequence
                  sb.append(currLine);
                  if ( newHeader.length() > 0 )
                  {
                     String[] array = newHeader.substring(1, newHeader.length()).split("\\s", 2);
                     if ( array.length >= 1 )
                        name = array[0];                  
                     if ( array.length >= 2 )
                        desc = array[1];
                     newHeader = "";
                  }
               }
            }
         }
      } catch (Exception e) { e.printStackTrace(); }

      return true;
   }

   public String getSequence()
   {
      return sequence;
   }

   public String getHeader()
   {
      return name + ' ' + desc;
   }

   public String getName()
   {
      return name;
   }

   public String getDesc()
   {
      return desc;
   }

   public int getNum()
   {
      return readNum;
   }

}
