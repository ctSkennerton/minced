import java.util.*;

public class CRISPR
{
   private Vector<Integer> repeats;
   private int repeatLength;

   public CRISPR()
   {
      repeats = new Vector<Integer>();
      repeatLength = 0;
   }

   public CRISPR(Vector<Integer> positions, int length)
   {
      repeats = positions;
      repeatLength = length;
   }

   public Vector repeats()
   {
      return repeats;
   }

   public int repeatLength()
   {
      return repeatLength;
   }

   public void setRepeats(Vector<Integer> _repeats)
   {
      repeats = _repeats;
   }

   public void setRepeatLength(int length)
   {
      repeatLength = length;
   }

   public int repeatSpacing(int pos1, int pos2)
   {
      return Math.abs(repeatAt(pos2) - repeatAt(pos1));
   }

   public void addRepeat(int val)
   {
      repeats.addElement(new Integer(val));
   }

   public void insertRepeatAt(int val, int pos)
   {
      repeats.insertElementAt(new Integer(val), pos);
   }

   public void setRepeatAt(int val, int pos)
   {
      repeats.setElementAt(new Integer(val), pos);
   }

   public void removeRepeat(int val)
   {
      repeats.removeElement(new Integer(val));
   }

   public int repeatAt(int i)
   {
      return ((Integer)repeats.elementAt(i)).intValue();
   }

   public int start()
   {
      return ((Integer)repeats.firstElement()).intValue();
   }

   public int end()
   {
      int lastRepeatBegin = ((Integer)repeats.lastElement()).intValue();
      return lastRepeatBegin + repeatLength - 1;
   }

   public int firstRepeat()
   {
      return ((Integer)repeats.elementAt(0)).intValue();
   }

   public int lastRepeat()
   {
      return ((Integer)repeats.lastElement()).intValue();
   }

   public int numRepeats()
   {
      return repeats.size();
   }

   public int numSpacers()
   {
      return numRepeats() - 1;
   }

   public String repeatStringAt(int i)
   {
      int currRepeatStartIndex = ((Integer)repeats.elementAt(i)).intValue();
      int currRepeatEndIndex = currRepeatStartIndex + repeatLength - 1;
      return DNASequence.seq.substring(currRepeatStartIndex, currRepeatEndIndex + 1);
   }

   public String spacerStringAt(int i)
   {
      int currRepeatEndIndex = ((Integer)repeats.elementAt(i)).intValue() + repeatLength - 1;
      int nextRepeatStartIndex = ((Integer)repeats.elementAt(i + 1)).intValue();
      int currSpacerStartIndex = currRepeatEndIndex + 1;
      int currSpacerEndIndex = nextRepeatStartIndex - 1;

      return DNASequence.seq.substring(currSpacerStartIndex, currSpacerEndIndex + 1);
   }

   public int averageSpacerLength()
   {
      int sum = 0;
      for (int i = 0; i < numSpacers(); i++)
      {
         sum = sum + spacerStringAt(i).length();
      }
      return sum/numSpacers();
   }

   public int averageRepeatLength()
   {
      int sum = 0;
      for (int i = 0; i < numRepeats(); i++)
      {
         sum = sum + repeatStringAt(i).length();
      }
      return sum/numRepeats();
   }

   public String toGff(String sequenceName, String parentName) {
	   String str = "";
	   for (int m = 0; m < numRepeats(); m++) {
		   int repeat_position = (repeatAt(m) + 1);
		   int repeat_end = repeat_position + this.repeatLength() - 1;
		   str += sequenceName + "\tminced:"+minced.VERSION+"\trepeat_unit\t" + repeat_position + "\t" + repeat_end + "\t1\t.\t.\tParent=" + parentName + ";ID=DR."+parentName+'.'+(m + 1) + "\n";
	   }
	   return str;
   }
   
   public String toString()
   {
      String str = "";

      String repeat, spacer, prevSpacer;
      repeat = spacer = prevSpacer = "";

      str += "POSITION\tREPEAT\t\t\t\tSPACER\n";

      str +="--------\t";

      for (int y = 0; y <  this.repeatLength(); y++)
         str +="-";
      str +="\t";

      for (int z = 0; z <  this.averageSpacerLength(); z++)
         str += "-";
      str +="\n";


      //add 1 to each position, to offset programming languagues that begin at 0 rather than 1
      for (int m = 0; m < numRepeats(); m++)
      {
         //repeat = getRepeat(m);
         str += (repeatAt(m) + 1) + "\t\t" + repeatStringAt(m) + "\t";

         // print spacer
         // because there are no spacers after the last repeat, we stop early (m < crisprIndexVector.size() - 1)
         if (m < numSpacers())
         {
            prevSpacer = spacer;
            spacer = spacerStringAt(m);
            str += spacer;

            str +="\t[ " + repeatStringAt(m).length() + ", " + spacerStringAt(m).length() + " ]";
            //str +="--[" + DNASequence.getSimilarity(repeatStringAt(m), spacerStringAt(m)) + "]";
            //str +="--[" + DNASequence.getSimilarity(spacer, prevSpacer) + "]";
            str += "\n";
         }
      }


      str +="\n--------\t";

      for (int x = 0; x < this.repeatLength(); x++)
         str += "-";
      str +="\t";

      for (int z = 0; z <  this.averageSpacerLength(); z++)
         str += "-";
      str +="\n";

      return str;
   }

}

