/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phylo.tree.reader.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * utilities for parsing data strings from file;
 * 
 * 
 * important notes about java String methods:
 *  1. endsWith() input string will not be treated as regular expression;
 *  2. split(), input string will be treated as regular expression;
 * 
 * @author tanxu
 */
public class StringUtils {
    public static final String EMPTY_STRING="";
    public static final String NEW_LINE_STRING = "\n"; //System.lineSeparator() == "\r\n";
    /**
     * split a data string s2 by delimiter s1 using the String.split method, with a modification for data string ends with delimiter;
     * 1. in the original split method, if the data string's ending substring is a n-copy of the delimiter, all the n empty strings will be removed;
     *  for example, datastring = "1abababab", delimiter = "ab", the split result will be [1]; datastring = "1abababab1", delimiter = "ab", the split result will be [1,"","","",1]
     * 2. this method will add the removed empty strings back to the resulted split array;
     * @param s1 regular expression
     * @param s2
     * @return 
     */
    public static List<String> splitS2ByS1AsExpected(String s1, String s2){
        List<String> ret;
        ret = CollectionUtils.copyArrayToList(s2.split(s1));
//        System.out.println("dataString:"+dataString);
//        System.out.println("delimiter:"+delimiter);
//        System.out.println(ret.size());
        
        //deal with ending delimiter regular expression copies
        int copyOfEndingDelimiters = copyOfS1AtEndOfS2(s1,s2);
//        System.out.println(copyOfEndingDelimiters);
        for(int i=0;i<copyOfEndingDelimiters;i++){
            ret.add(EMPTY_STRING);
        }
        
        //deal with begnining empty string if the full data string is composed of delimiters
        if(stringS2IsFullyDecomposableToS1(s1,s2)){
            ret.add(EMPTY_STRING);
        }
        
//        System.out.println("done");
        return ret;
    }
    
    
    
    
    
    public static List<String> splitStringByListOfSubstringLenghts(List<Integer> listOfLengths, String dataString){
        List<String> ret = new ArrayList<>();
        
        
        
        //todo
        return ret;
    }
    
    
    /**
     * number of duplicates of s1 at the end of s2;
     * for example s2 = "1111111", s1 = 1, copy number  = 7;
     * s2 = "sdfaabcabfabcab2", s1 = "ab.", copy number = 4
     * @param s1regex regular expression
     * @param s2 data string
     * @return 
     */
    public static int copyOfS1AtEndOfS2(String s1regex, String s2){
        int copy = 0;
        Pattern p = Pattern.compile("(.*)(".concat(s1regex).concat(")"));//group 1 is the leading string, group 2 is the ending string
//        System.out.println(p);
        Pattern p1 = Pattern.compile("(.*)(\\n)");
//        System.out.println(p1);
        
        Matcher m = p.matcher(s2);
        while(m.matches()){
            copy++;
            m = p.matcher(m.group(1));
        }
        
        
        
        return copy;
    }
    
    /**
     * whether or not string s2 is fully a copy of duplicates of regular expression s1regex;
     * that is to say, s2 is from beginning to end, can be decomposed to individual s1regex regular expressions;
     * @param s1regex
     * @param s2
     * @return 
     */
    public static boolean stringS2IsFullyDecomposableToS1(String s1regex, String s2){
        Pattern p = Pattern.compile("(.*)(".concat(s1regex).concat(")"));//group 1 is the leading string, group 2 is the ending string
        Matcher m = p.matcher(s2);
        
        String remainingString=s2;
        while(m.matches()){
            remainingString = m.group(1);
            m = p.matcher(remainingString);
        }
        
        return remainingString.isEmpty();
    }
    
    
    
    
    
    /**
     * split data string by s1 that is not following a s2;
     * the dataString might contain substrings s1 o s2s1, but only split the data string by s1 that is not following a s2;
     * be cautious if s1 and/or s2 are regular expression string with quantifiers;
     * for example, 
     *  if dataString = "1111=231231\\=321323\\= "; s1 = "=", s2 = "\\"; the resulting split should be ("1111", "231231\\=321323\\= ");
     * @param s1
     * @param s2
     * @param dataString 
     * @return  
     */
    public static List<String> splitByS1NotFollowingS2(String s1, String s2, String dataString){
        
        String s21 = s2.concat(s1);
//        System.out.println(s21);
        
        if(dataString.contains(s21)){
//            System.out.println("containing s21");
            List<String> correctSplits = new ArrayList<>();
            List<String> splits0 = splitS2ByS1AsExpected(s1, dataString);
            String recoveredString = "";
            for(String s:splits0){
//                System.out.println("split:"+s);
                recoveredString = recoveredString.concat(s);
                String testString = recoveredString.concat(s1);
//                System.out.println("testString:"+testString);
                if(testString.endsWith(s21)){//need to recover
//                    System.out.println("ends with s21");
                    recoveredString = recoveredString.concat(s1);
                }else{//correct split
//                    System.out.println("correct split:"+recoveredString);
                    correctSplits.add(recoveredString);
                    recoveredString = "";
                }
            }
            
            
            return correctSplits;
            
        }else{
             return splitS2ByS1AsExpected(s1,dataString);
        }
        
    }
    
    
    /**
     * 
     * @param delimiterList
     * @param dataString
     * @return 
     */
    public static List<String> splitByAListOfDelimiters(List<String> delimiterList, String dataString){
        
        String line = "This order was placed for QT3000! OK?";
      String pattern = "(.*)(\\d+)(.*)";

      // Create a Pattern object
      Pattern r = Pattern.compile(pattern);

      // Now create matcher object.
      Matcher m = r.matcher(line);
      
      if (m.find( )) {
          System.out.println(m.groupCount());
         System.out.println("Found value: " + m.group(0) );
         System.out.println("Found value: " + m.group(1) );
         System.out.println("Found value: " + m.group(2) );
          System.out.println(m.group(3));
      } else {
         System.out.println("NO MATCH");
      }
        return null;
        
        //todo
    }
    
    /**
     * split the data string by the first n delimiters; all the remaining string that containing delimiter will be put together in the last split component the same as the original string
     * @param delimiter
     * @param dataString
     * @param firstN
     * @return 
     */
    public static List<String> splitByFirstNDelimiters(String delimiter, String dataString, int firstN){
        
        List<String> splits = splitS2ByS1AsExpected(delimiter, dataString);
//        System.out.println("=====splits:"+dataString);
//        printStringList(splits);
//        System.out.println("=====");
        
        if(firstN+1>=splits.size()){
//            System.out.println("direct return");
            return splits;
        }else{//there are more than firstN delimiters in the data string
            List<String> ret = new ArrayList<>();
            String regex = "";
            for(int i=0;i<firstN;i++){
                regex = regex.concat(splits.get(i)).concat(delimiter);
                ret.add(splits.get(i));
            }
//            System.out.println("regex:"+regex);
            Pattern p = Pattern.compile(regex);
            Matcher m = p.matcher(dataString);
            
            if(m.find()){
//                System.out.println("found");
//                System.out.println(m.end(0));
//                System.out.println(dataString.substring(m.end(0)));
                ret.add(dataString.substring(m.end(0)));
//                System.out.println(m.end(1));
            }
//            StringUtils.printStringList(ret);
            return ret;
            
        }
        
        
    }
    
    /**
     * s2 is a string that matches regular expression string s1;
     * for example: s1 is "\\s+", s2 is "\t\n ", it will return true;
     * @param s1
     * @param s2
     * @return 
     */
    public static boolean stringS2IsCoveredByS1AsRegex(String s1, String s2){
        Pattern p = Pattern.compile(s1);
        return p.matcher(s2).matches();
    }
    
    
    
    /**
     * data string s2 contains a substring that matches regular expression string s1
     * @param s1
     * @param s2
     * @return 
     */
    public static boolean stringS2ContainsS1AsRegex(String s1, String s2){
        return Pattern.compile(s1).matcher(s2).find();
    }
    
    
    /**
     * check when an input string is processed with some string method with the given s1 as parameter, whether the given string s2 will interfere with the processing or not;
     * return true if s1 can be successfully processed without s2 being involved, false otherwise;
     * <p>
     * for example, let s1 be record delimiter, s2 be mandatory attribute delimiter, when split a data string with s1 into record strings, s2 will not be split;
     * </p><p>
     * that is to say, there should not be any substring in s2 that is recognized as s1 by the specific string method;</p>
     * 
     * <p>
     * precedence list:
     * metadata identifier > (if has metadata)
     * segment delimiter (if segmental structure)
     * record delimiter >
     * attribute delimiter > 
     * mandatory attribute null value string > 
     * composite tag Attribute Delimiter (if composite tag attribute exists)
     * tag attribute component delimiter (if composite tag and/or tailing tag attribute exist)
     * tag attribute component multi-value delimiter (if allowing multi-value)
     * 
     * ...
     * 
     * </p>
     * @param s1
     * @param s2
     * @param toUseSplit if true, the related string process method is split(String regex), the s1 and s2 are used as regular expression;
     * @param toUseContains if true, the related string process method is contains(String str), the s1 and s2 are used as plain string;
     * @return true if s1 can be successfully processed without s2 being involved, false otherwise;
     */
    public static boolean s1IsHigherPriorityThanS2(String s1, boolean s1IsPlainString, String s2, boolean s2IsPlainString, boolean toUseSplit, boolean toUseContains){
        //
        
//        if(toUseSplit&&toUseContains){
//            throw new IllegalArgumentException("toUseSplit and toUseContains cannot both be true");
//        }
//        if(!toUseSplit && !toUseContains){
//            throw new IllegalArgumentException("toUseSplit and toUseContains cannot both be false");
//        }
//        
//        if(toUseSplit){//regular expression
//           return stringS2ContainsS1AsRegex(s1, s2);
//        }
//        
//        
//        if(toUseContains){
//            return s2.contains(s1);
//        }
//        
        
        return true;
        
    }
    
    
    
    /**
     * test if the given string contains only white spaces;
     * @param testString
     * @return 
     */
    public static boolean isWhiteSpaces(String testString){
//        testString.
        Pattern p = Pattern.compile("\\s+");
        Matcher m = p.matcher(testString);
        return m.find();
    }
    
    
    
    
    
    
    /**
     * test if the given string is empty
     * @param testString
     * @return 
     */
    public static boolean isEmptyString(String testString){
        return testString.isEmpty();
    }
    
    
    /**
     * a notable java string split() method's specification:
     * 1. if the data string is ended with the delimiter, the last delimiter will be removed from data string before the data string is split;
     * thus, the resulting array will not contain a empty string at the end;
     * 2. but if the data string is started with delimiter, the first element of the resulting array will be an empty string;
     */
    public static void splitTest(){
//        String s1 = "0;1;2;3;4; ";
//        String[] split1 = s1.split(";");
//        printStringArray(split1);
        
        
//        String s2 = "0;1;2;3;4;";
//        String[] split2 = s2.split(";");
//        printStringArray(split2);
        
        
        //if data string ends with delimiter, no split element will be added in the array after the last delimiter; 
        //however, if two delimiters are adjacent to each other, an empty string will be added as a split elements in the array
//        String s3 = "0;1;2;3;;4;";
//        String[] split3 = s3.split(";");
//        printStringArray(split3);
        
        
        String s4 = " ;0;1;2;3;;4";
        String[] split4 = s4.split(";");
        printStringArray(split4);
        
        
        String s5 = ";0;1;2;3;;4";
        String[] split5 = s5.split(";");
        printStringArray(split5);
    }
    
    
    public static void printStringArray(String[] ss){
        for(int i=0;i<ss.length;i++){
            System.out.println(i+":"+ss[i]);
        }
    }
    public static void printStringList(List<String> ss){
        for(int i=0;i<ss.size();i++){
            System.out.println(i+":"+ss.get(i));
        }
    }
    
    
    public static void catenateTest(){
        String cat = "\t+";
        String s1 = "1111";
        String s2 = "2222";
        String s3 = "           ";
        System.out.println(s3.equals(cat));
        System.out.println(s3.contains(cat));
        System.out.println(s1.concat(cat).concat(s2));
    }
    
    public static void main(String[] args){
//        System.out.println("hello");
        
//        String s = "1=2=3=4=";
//        String s1 = "=1=2=3=4=";
//        String s2 = "=1==2=3=4=";
//        String s3 = "1111=231231\\=321323\\= ";
        
//        String[] splits = s.split("=");
//        System.out.println(s.split("=").length);
//        System.out.println(s1.split("=").length);
//        System.out.println(s2.split("=").length);
//        printStringArray(s3.split("="));
        
        
//        splitByS1NotFollowingS2("=","\\",s3);
//        splitTest();
        
//        System.out.println(isWhiteSpaces("\t\n\r\f"));
//        System.out.println(isEmptyString("\t"));


//            splitByAListOfDelimiters(null,null);
//            catenateTest();
//            String delimiter = "\t+";
//            String data = "1kdfa\tjiijfds\t\t\tkjdfa\t\t\t\t\tdfadsfai\t ddd";
//            
//            String delimiter2 = "\\s+";
//            String data2 = "1kdfa\t   jiijfds\t\t\t\nkjdfa\t\t\t\t    \tdfadsfai            \t ddd";
//            printStringList(splitByFirstNDelimiters(delimiter2, data2, 2));

//            String s = "dfas;jdfsaj,(dfafadfadsfadfadff),(dssss   ";
//            
//            printStringArray(s.split("\\,\\("));
            
//            String s2 = "sdfaabcabfabcab2";
//            String s1 = "ab.";
//            String s = "1sssssdsssssssss";
//            String delimiter = "s+";
//            System.out.println("hello");
//            printStringArray(s2.split(s1));
//            System.out.println(s2.endsWith(s1));
            
//            List<String> ret = Arrays.asList(s.split(delimiter));
//            ret.add("");

//            String s2 = "/";
//            String s1 = "/";
//            System.out.println(copyOfS1AtS2Ending(s1,s2));
//            printStringList(splitAsExpected(s1,s2));
            
//            System.out.println(S2IsFullyDecomposableToS1(s1, s2));
                    
        String s1 = "\\s+";
        String s2 = "\t     \n          ";
        System.out.println(stringS2IsCoveredByS1AsRegex(s1, s2));
    }
    
    
}
