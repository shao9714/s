import java.io.*;
import java.lang.String;
import java.lang.Object;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.regex.*;
import org.apache.commons.math3.distribution.*;

public class squant {

    static String inputFile = new String(); 
    static String outputFile = new String();

    static int transcriptCount = 0;
    static Map<String, Integer> transcripts = new LinkedHashMap<>();
    static Map<Integer, String[][]> reads = new LinkedHashMap<>();
    
    static double mean=200, sd=25;
    static NormalDistribution dist = new NormalDistribution(mean, sd);

    static Map<String, Double> trans_length = new LinkedHashMap<>();
    static Map<String, Double> trans_count = new LinkedHashMap<>();
    static Map<Integer, Double[][]> reads_p = new LinkedHashMap<>();
    static Map<String, Double> trans_p1 = new LinkedHashMap<>();
    static Map<String, Double> trans_p1_prev = new LinkedHashMap<>();
    
    public static void main(String args[]) {
        int input=-1, output=-1;
        
        if (args.length < 4) {
            System.out.println("Not enough input variables!");
        } else if (args.length > 5) {
            System.out.println("Too many input variables!");
        } else if (args.length == 5 && !args[4].equals("--eqc")) {
            System.out.println("Invalid input variable(s)!");
        } else {
            if (args[0].equals("--in") && args[2].equals("--out")) {
                input = 1;
                output = 3;
            } else if (args[0].equals("--out") && args[2].equals("--in")) {
                input = 3;
                output = 1;
            } else {
                System.out.println("Invalid input variable(s)!");
            }
            inputFile = args[input];
            outputFile = args[output];
            readFile();
            System.out.println("Done reading file.");
            if (args.length == 5 && args[4].equals("--eqc")) {
                System.out.println("eqc() not implemented!");
            } else {
                fullEM();
                System.out.println("Done running fullEM");
                writeFile();
                System.out.println("Done writing file");
            }
        }
    }

    public static void fullEM() {
        boolean notConverged = true;
        int iter = 0;
        
        // Storing all effective lengths
        storeEffL();
        System.out.println("Done storing effective lengths");
        
        // Storing all constant p2, p3, p4 for all reads
        setUpReadsP();
        System.out.println("Done setting up reads");
        System.out.println(reads.size());
        
        // EM Algorithm starts here
        // First E step
    	initialEStep();
    
        while (notConverged) {
        	iter++;
        	System.out.println(iter + "th Round");
        	
            // M step
            for (Map.Entry<String, Double> entry : trans_count.entrySet()) {
            	trans_p1.put(entry.getKey(), entry.getValue()/reads.size());
            }
            
            if (converged()) {
            	System.out.println("Converged");
                notConverged = false;
                break;
            }
           
            for (Map.Entry<String, Double> entry : trans_count.entrySet()) {
            	trans_count.put(entry.getKey(), 0.0);
            }
            
            for (Map.Entry<String, Double> entry : trans_p1.entrySet()) {
            	trans_p1_prev.put(entry.getKey(), entry.getValue());
            }
            
            // Repeat E Step
            for (Map.Entry<Integer, Double[][]> entry : reads_p.entrySet()) {
            	Integer index = entry.getKey();
                Double[][] p = entry.getValue();
                Double[] sum = new Double[p.length];
                double norm=0;
                String[][] tNameA = reads.get(index);
                
                for (int i=0; i<p.length; i++) {
                	String tName = tNameA[i][0];
                	sum[i] = trans_p1.get(tName)*p[i][1]*p[i][2]*p[i][3];
                	norm+=sum[i];
                }
                for (int i=0; i<p.length; i++) {
                	sum[i]/=norm;
                }
                
                String[][] updateTCount = reads.get(index);
                
                for (int i=0; i<updateTCount.length; i++) {
                	String name = updateTCount[i][0];
                	double currCount = trans_count.get(name);
                	currCount += sum[i];
                	trans_count.put(name, currCount);
                }
                
            }
        }
    }

    public static boolean converged() {
    	for (Map.Entry<String, Double> entry : trans_p1.entrySet()) {
    		String name = entry.getKey();
    		Double val = entry.getValue();
    		
    		if (!equalUpTo(val, trans_p1_prev.get(name))) {
    			return false;
    		}
    	}
    	
    	return true;
    }
    
    public static boolean equalUpTo(Double a, Double b) {
    	if (a < 0.00001 && b < 0.0001) {
    		return true;
    	}
    	
    	String aString = a.toString();
    	String bString = b.toString();
    	int length = aString.length() > bString.length() ? bString.length() : aString.length();
    	int count = 0;
    	
    	int i, j;
    	
    	for (i=0; i<aString.length(); i++) {
    		if (aString.charAt(i)=='.') {
    			i++;
    			break;
    		}
    	}
    	
    	for (j=0; j<bString.length(); j++) {
    		if (bString.charAt(j)=='.') {
    			j++;
    			break;
    		}
    	}
    	
    	for (int k=0; k<length; k++) {
    		if (count==2) {
    			return true;
    		}
    		if (k>3 || i==aString.length() || j==bString.length()) {
    			return false;
    		}
    		if (aString.charAt(i)=='E' || bString.charAt(j)=='E') {
    			return false;
    		}
    		
    		if (aString.charAt(i)==bString.charAt(j)) {
    			count++;
    		}
    		i++;
    		j++;
    	}
    
    	return false;
    }
    
    public static void initialEStep() {
    	 Double[] sum;
    	 int iter = 0;
    	 
    	 for (Map.Entry<Integer, Double[][]> entry : reads_p.entrySet()) {
         	 Integer index = entry.getKey();
             Double[][] p = entry.getValue();
             sum = new Double[p.length];
             double norm=0.0;
             iter++;
             
             if (iter % 2000==0) {
            	 System.out.println("Initial run at " + iter);
             }
             
             for (int i=0; i<p.length; i++) {
             	sum[i] = p[i][0]*p[i][1]*p[i][2]*p[i][3];
             	norm+=sum[i];
             }
             
             for (int i=0; i<p.length; i++) {
             	sum[i]/=norm;
             }
             
             for (int i=0; i<p.length; i++) {
            	 String[][] updateTCount = reads.get(index);
            	 String name = updateTCount[i][0];
            	 double currCount = trans_count.get(name);
            	 currCount += sum[i];
            	 trans_count.put(name, currCount);
             }
         }
         
    	 
    	 for (Map.Entry<String, Double> e : trans_count.entrySet()) {
    		 double d = 1.0/reads.size();
    		 trans_p1_prev.put(e.getKey(), d);
         }
    }
    
    public static void storeEffL() {
    	for (Map.Entry<String, Integer> entry : transcripts.entrySet()) {
        	String name = entry.getKey();
        	int length = entry.getValue();
        	double d = effectiveLength(length);
        	trans_length.put(name, d);
        }
    }
    
    public static void setUpReadsP() {
    	double p1=1.0/reads.size(), p2, p3, p4;
    	int iter = 0;
    	System.out.println("Setting up reads.");
    	for (Map.Entry<Integer, String[][]> entry : reads.entrySet()) {
        	Integer index = entry.getKey();
        	String[][] alignment = entry.getValue();
        	Double[][] val= new Double[alignment.length][4];
            String name, ori; 
            int pos;
           
            
            iter++;
        	if (iter%2000==0) {
        		System.out.println("Set-up at read " + iter);
        	}
            
            for (int i=0; i<alignment.length; i++) {
            	
                name = alignment[i][0];
                ori = alignment[i][1];
                pos = Integer.valueOf(alignment[i][2]);
                p4 = Double.valueOf(alignment[i][3]);
                
                // Extracting effLength
                int length = transcripts.get(name);
              	double d =  trans_length.get(name);
                double effLength = d;
                //

                p2 = 1/effLength;
                if (ori.equals("f")) {
                	p3 = dist.cumulativeProbability(length-pos);
                } else {
                	p3 = dist.cumulativeProbability(pos+100);
                }
                
                val[i][0] = p1;
                val[i][1] = p2;
                val[i][2] = p3;
                val[i][3] = p4;
            }
            reads_p.put(index, val);
        }
    }
    
    public static double effectiveLength(int length) {

        if (length > 1000) {
            return length-200;
        } else {
            double sum=0, condMean=0;
            double[] p = new double[length+1];
            for (int i=0; i<=length; i++) {
                p[i]=dist.density(i);
                sum+=p[i];
            }
            for (int i=0; i<p.length; i++) {
                p[i]/=sum;
            }
            for (int i=0; i<p.length; i++) {
                condMean += i*p[i];
            }

            return length-condMean;
        }
    }

    public static void writeFile() {
        BufferedWriter bw;
        FileWriter fw;

        try {
            File file = new File(outputFile);

            if (!file.exists()) {
                file.createNewFile();
            }

            fw = new FileWriter(file);
            bw = new BufferedWriter(fw);
        
            for (Map.Entry<String, Double> entry : trans_count.entrySet()) {
            	bw.write(entry.getKey() + "\t" + (double)Math.round(trans_length.get(entry.getKey())) + 
            			"\t" + (double)Math.round(entry.getValue()) + "\n");
            }

            bw.close();
            fw.close();
        } catch (IOException e) {
            System.out.println("File writing error " + outputFile + " !");
            e.printStackTrace();
        }
    }

    public static void readFile() {
        FileReader fr; 
        BufferedReader br;
        String str;
        String[][] list;
        int count=0, alnCount=0;
        Pattern pattern = Pattern.compile("^(.*)[\t](.*)$");
        Matcher matcher;

        try {
            fr = new FileReader(inputFile);
            br = new BufferedReader(fr);
            
            transcriptCount = Integer.valueOf(br.readLine());
            
            while ((str = br.readLine()) != null && !str.matches("^[0-9]*$")) {
                matcher = pattern.matcher(str);
                if (matcher.find()) {
                    transcripts.put(matcher.group(1), Integer.valueOf(matcher.group(2)));
                    trans_length.put(matcher.group(1), 0.0);
                    trans_count.put(matcher.group(1), 0.0);
                } else {
                    System.out.println("Error in readFile(). Wrong format in transcript records.");
                    br.close();
                    fr.close();
                    return;
                }
            }

            int c = 0;
            while (str != null) {
                list = new String[Integer.valueOf(str)][4];
                while ((str = br.readLine()) != null && !str.matches("^[0-9]*$")) {
                    list[alnCount++] = str.split("[\t]");
                }
                reads.put(count++, list);
                alnCount=0;
                if (str == null) {
                    break;
                }
                c++;
                if (c%2000 ==0 ) {
                	System.out.println("Reading reads in progress " + c);
                }
            }

            br.close();
            fr.close();
        } catch (Exception e) {
            System.out.println(e.getClass());
            System.out.println("File reading error " + inputFile + "!");
        }
    }

}