import java.io.*;

public class minced
{
   public static void main(String[] args)
   {
      //default values
      int screenDisplay      = 1;
      int minNumRepeats      = 3;
      int minRepeatLength    = 19;
      int maxRepeatLength    = 38;
      int minSpacerLength    = 19;
      int maxSpacerLength    = 48;
      int searchWindowLength = 8;
      
      int numOptions = 0;

      if (args.length == 0)
      {
         printUsage();
         //System.out.println("No arguments specified.");
         //System.out.println("Try 'minced --help' for more information.");
         System.exit(1);
      }

      if (args[0].equals("--help"))
      {
         printUsage();
         System.exit(1);
      }
      
      if (args[0].equals("--version"))
      {
         printVersion();
         System.exit(1);
      }
            
      int i = 0;
      while (args[i].charAt(0) == '-')
      {
         try
         {
            if (args[i].endsWith("minNR"))
            {
               i++;
               if (i >= args.length)
               {
                  System.out.println("Missing option value for " + args[i - 1] + ".");
                  System.out.println("Try 'minced --help' for more information.");
                  System.exit(1);
               }
               minNumRepeats = Integer.parseInt(args[i]);
               numOptions++;
            }
            else if (args[i].endsWith("minRL"))
            {
               i++;
               if (i >= args.length)
               {
                  System.out.println("Missing option value for " + args[i - 1] + ".");
                  System.out.println("Try 'minced --help' for more information.");
                  System.exit(1);
               }
               minRepeatLength = Integer.parseInt(args[i]);
               numOptions++;
            }
            else if (args[i].endsWith("maxRL"))
            {
               i++;
               if (i >= args.length)
               {
                  System.out.println("Missing option value for " + args[i - 1] + ".");
                  System.out.println("Try 'minced --help' for more information.");
                  System.exit(1);
               }
               maxRepeatLength = Integer.parseInt(args[i]);
               numOptions++;
            }
            else if (args[i].endsWith("minSL"))
            {
               i++;
               if (i >= args.length)
               {
                  System.out.println("Missing option value for " + args[i - 1] + ".");
                  System.out.println("Try 'minced --help' for more information.");
                  System.exit(1);
               }
               minSpacerLength = Integer.parseInt(args[i]);
               numOptions++;
            }
            else if (args[i].endsWith("maxSL"))
            {
               i++;
               if (i >= args.length)
               {
                  System.out.println("Missing option value for " + args[i - 1] + ".");
                  System.out.println("Try 'minced --help' for more information.");
                  System.exit(1);
               }
               maxSpacerLength = Integer.parseInt(args[i]);
               numOptions++;
            }               
            else if (args[i].endsWith("searchWL"))
            {
               i++;
               if (i >= args.length)
               {
                  System.out.println("Missing option value for " + args[i - 1] + ".");
                  System.out.println("Try 'minced --help' for more information.");
                  System.exit(1);
               }
               searchWindowLength = Integer.parseInt(args[i]);
               numOptions++;
               if ( (searchWindowLength < 6) || (searchWindowLength > 9) )
               {
                  searchWindowLength = 8;
                  System.out.println("--Reseting search window length to " + searchWindowLength + "--");
               }
            }
            /*
            else if (args[i].endsWith("screen"))
            {
               i++;
               if (i >= args.length)
               {
                  System.out.println("Missing option value for " + args[i - 1] + ".");
                  System.out.println("Try 'minced --help' for more information.");
                  System.exit(1);
               }
               screenDisplay = Integer.parseInt(args[i]);
               numOptions++;
               if ( screenDisplay != 1 )
                  screenDisplay = 0;
            }
            */
            else
            {
               System.out.println("Invalid option:  " + args[i]);
               System.out.println("Try 'minced --help' for more information.");
               System.exit(1);
            }
            i++;  //--------//
            if (i >= args.length)
            {
               System.out.println("No input file specified.");
               System.out.println("Try 'minced --help' for more information.");
               System.exit(1);
            }
         }  //end try
         catch (NumberFormatException e) // exception will be caught if user does not specify an integer after each option
         {
            System.out.println("Invalid argument type.  Please, specify an option followed by an integer.  For example:  -minNR 3");
            System.out.println("Try 'minced --help' for more information.");
            System.exit(1);
         }
      }  // end while
      
      
      // Last options should be an input file and optional output file
      String inputFileName = "";
      String outputFileName = "";
      boolean outputFileSpecified = false;
      int numArgsRemaining = args.length - (2 * numOptions);
      
      if (numArgsRemaining == 1)
         inputFileName = args[i];
      else if (numArgsRemaining == 2)
      {
         inputFileName = args[i];
         outputFileSpecified = true;
         outputFileName = args[i + 1];
         screenDisplay = 0;
      }
      else
      {
         System.out.println("Improper usage.");
         System.out.println("Try 'minced --help' for more information.");
         System.exit(1);
      }
      
      
      File inputFile = new File(inputFileName);
      if (!inputFile.exists())
      {
         System.out.println("Input file " + inputFile.getPath() + " does not exist.");
         System.exit(1);
      }

      if (inputFile.isDirectory())
      {
         System.out.println("You have entered a directory name. An input file name is required: " + inputFile.getPath());
         System.exit(1);
      }


      
      File outputFile;
      if (outputFileSpecified)
      {
         outputFile = new File(outputFileName);
         
         if (outputFile.isDirectory())
         {
            System.out.println("You have entered an existing directory name. An output file name is required: " + outputFile.getAbsolutePath());
            System.exit(1);
         }
         
         if (!outputFile.getAbsoluteFile().getParentFile().isDirectory())
         {
            System.out.println("You did not enter a valid file output path: " + outputFile.getAbsolutePath());
            System.exit(1);
         }
      }


      CRISPRFinder client = new CRISPRFinder(inputFileName, outputFileName, screenDisplay, minNumRepeats, minRepeatLength, maxRepeatLength, minSpacerLength, maxSpacerLength, searchWindowLength);
      client.goCRISPRFinder();

   }
   
   public static void printUsage()
   {  
      System.out.println("MinCED, a program to find CRISPRs in shotgun DNA sequences or full genomes");
      System.out.println();
      System.out.println("Usage:    java minced [options] inputMultiFastaFile [outputFile]");
      System.out.println();
      System.out.println("Options:  -searchWL  Length of search window used to discover CRISPRs (range: 6-9). Default: 8");
      System.out.println("          -minNR     Minimum number of repeats a CRISPR must contain. Default: 3");
      System.out.println("          -minRL     Minimum length of the CRISPR repeats. Default: 19");
      System.out.println("          -maxRL     Maximum length of the CRISPR repeats. Default: 38");
      System.out.println("          -minSL     Minimum length of the CRISPR spacers. Default: 19");
      System.out.println("          -maxSL     Maximum length of the CRISPR spacers. Default: 48");
      //System.out.println("          -screen    Print results to the screen, instead of a file; (range: 0-1); default 0");
      System.out.println();
      System.out.println("Examples: java minced ecoli.fna");
      System.out.println("          java minced -minNR 2 metagenome.fna");
      System.out.println("          java minced -minNR 2 metagenome.fna metagenome.crisprs");
      System.out.println();
   }

   public static void printVersion()
   {  
      System.out.println("MinCED - Mining CRISPRs in Environmental Datasets (version 0.1)");
      System.out.println("Copyright 2011 Florent ANGLY <florent.angly@gmail.com>");
      System.out.println("Distributed under the GNU General Public License version 3");
      System.out.println();
   }

}