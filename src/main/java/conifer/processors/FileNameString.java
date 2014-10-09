package conifer.processors;

public class FileNameString
{
  /*
   * This code serves to extract the number of sites in the file names
   * this number of sites is the number of sites in the alignment of the input data
   */
  private String fileName;
  public FileNameString()
  {
    this.fileName=null;
  }
  public FileNameString(String fileName)
  {
    this.fileName = fileName;
  }
  
  public String extractDigits(String src) {
    StringBuilder builder = new StringBuilder();
  for (int i = 0; i < src.length(); i++) {
      char c = src.charAt(i);
      if (Character.isDigit(c)) {
        builder.append(c);
      }
  }
  return builder.toString();
 }
  

  public String subStringBetween(String sentence, String before, String after) {

        int startSub = this.subStringStartIndex(sentence, before);
        int stopSub = this.subStringEndIndex(sentence, after);

        String newWord = sentence.substring(startSub, stopSub);
        return newWord;
    }

    public int subStringStartIndex(String sentence, String delimiterBeforeWord) {

       int startIndex = sentence.lastIndexOf(delimiterBeforeWord);
       int result = startIndex+delimiterBeforeWord.length(); 
       return (result);
    }
    
 
    public int subStringEndIndex(String sentence, String delimiterAfterWord) {

        int startIndex = 0;
        String newWord = "";
        int x = 0;

        for (int i = 0; i < sentence.length(); i++) {
            newWord = "";

            if (sentence.charAt(i) == delimiterAfterWord.charAt(0)) {
                startIndex = i;
                for (int j = 0; j < delimiterAfterWord.length(); j++) {
                    try {
                        if (sentence.charAt(startIndex) == delimiterAfterWord.charAt(j)) {
                            newWord = newWord + sentence.charAt(startIndex);
                        }
                        startIndex++;
                    } catch (Exception e) {
                    }

                }
                if (newWord.equals(delimiterAfterWord)) {
                    x = startIndex;
                    x = x - delimiterAfterWord.length();
                }
            }
        }
        return x;
    }
  
  public static void main(String [] args){
    String fileName = "ABCnumberOfSites100Seeds100";
    FileNameString fileNameString = new FileNameString(fileName);
    System.out.print(fileNameString.subStringBetween(fileName, "numberOfSites", "Seeds"));
    String numInFileName = fileNameString.extractDigits(fileName);
    String segments[]= fileName.split("numOfSites");
  }  

}


