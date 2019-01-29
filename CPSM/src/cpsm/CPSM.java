/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cpsm;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.docopt.Docopt;

/**
 *
 * @author minxian
 */
public class CPSM {
    
    private static final LinkedList<String[]> ContentList = new LinkedList<>();
    private static double [] P1_Array = null;
    private static double [] P2_Array = null;
    private static int[] pv_index = null;       //Column index for p values.
    private static int POS_index = -1;          //Column index for physical position.
    private static int[]  POS_Array = null;
    private static int[][] WINDON_BOUNDARY = null;
    private static int half_windown_size = 30 * 1000;
    private static int Permute_P1 = -1;
    private static int Permute_P2 = -1;
    private static double[] results = null;
    private static final DecimalFormat formator = new DecimalFormat(".####");
    
    
    private static final String DOC =
                "CPSM method for checking the shared driver snps between a pair of GWAS.\n"
                + "\n"
                + "Usage:\n"
                + "  CPSM -p int,int -b int -m int [--p1 int | --p2 int] [--half_window_size int] [-t cpus]\n"
                + "  CPSM (-h | --help)\n"
                + "  CPSM --version\n"
                + "\n"
                + "---------------------------\n"
                + "Read data from stdin, output results to stdout.\n"
                + "Each non-title row is a result, multiple rows for permutation test.\n"
                + "---------------------------\n"
                + "\n"
                + "Options:\n"
                + "  -p int,int    Column index for pvalues, two ints, index starts from 1.!\n"
                + "  -b int        Column index for physical position.\n"
                + "  -m int        Column index for mark id, output as title.\n"
                + "  --half_window_size int    Half window size in kb, default 30. \n"
                + "  --p1 int      Permute the first summary int times, close the output for original test.\n"
                + "  --p2 int      Permute the second summary int times, close the output for original test.\n"
                + "  -t cpus       Number of cpus for computing.\n"
                + "  -h --help     Show this screen.\n"
                + "  --version     Show version.\n"
                + "\n";
    
    
    private static void setWindowBoundary(){
        WINDON_BOUNDARY = new int[POS_Array.length][2];
        int l = 0;
        int r = 0;
        for (int f = 0; f < P1_Array.length; f++) {
            while (POS_Array[f] - POS_Array[l]>half_windown_size) {
                l++;
            }
            while (r < POS_Array.length && POS_Array[r] - POS_Array[f]<half_windown_size) {
                r++;
            }
            
            WINDON_BOUNDARY[f][0] = l;
            WINDON_BOUNDARY[f][1] = r;
        }
    }
    
    /**
     * Compute the CPSM sore at a single position.
     * @param focal
     * @param left
     * @param right
     * @return 
     */
    private static double getSingleScore(int focal, int left, int right){
        if (left == right) {
            return 0.0;
        }
        
        int win_len = POS_Array[right] - POS_Array[left];
        double[] w = IntStream.range(left, right+1) //include the right end.
                .mapToDouble(i->{
                    return 1 - Math.abs(POS_Array[focal]-POS_Array[i])*1.0/win_len;
                })
                .toArray();
       
        double sum_w = Arrays.stream(w).sum();
        double[] w_normal = Arrays.stream(w).map(i->i/sum_w).toArray();
        double top = IntStream.range(left,right+1)
                .mapToDouble(i -> {
                    return P1_Array[i] * P2_Array[i] * w_normal[i];
                })
                .sum();
        double bottom = 1- Arrays.stream(w_normal).map(i->i*i).sum();
                
        if (bottom == 0) {
            return -1000;
        }else{
            return top/bottom;
        }
    }
    
    private static void outputAllscores(){
        IntStream.range(0, WINDON_BOUNDARY.length)
                .parallel()
                .forEach(i -> {
                    results[i] = getSingleScore(i, WINDON_BOUNDARY[i][0], WINDON_BOUNDARY[i][1]);
                });
        
        String s = Arrays.stream(results)
                .mapToObj(i -> formator.format(i))
                .collect(Collectors.joining("\t"));
        System.out.println(s);
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        Map<String, Object> opts =
                 new Docopt(DOC).withVersion("1.8").parse(args);
        
        if(opts.get("-p") != null){
               String[] temp = ((String) opts.get("-r")).split(",");
               pv_index = Arrays.stream(temp)
                       .mapToInt(s->Integer.parseInt(s)-1)
                       .toArray();
        }
        
        if(opts.get("-b") != null){
               POS_index = Integer.parseInt((String) opts.get("-b"));
        }
        if(opts.get("--half_window_size") != null){
               half_windown_size = Integer.parseInt((String) opts.get("--half_window_size")) * 1000;
        }
        if(opts.get("--p1") != null){
               Permute_P1 = Integer.parseInt((String) opts.get("--p1"));
        }
        if(opts.get("--p2") != null){
               Permute_P2 = Integer.parseInt((String) opts.get("--p2"));
        }
        if(opts.get("-t") != null){   
            System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", (String) opts.get("--p2"));
        }
        int MARKID_COL = Integer.parseInt((String) opts.get("-m"))-1;
        
        
        BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
            in.lines()
              .filter(s -> s.length() > 0)
              .sequential() //*** Very Important **** Here must to be sequential, because we will update index.
              .forEach(s ->{
                  String[] ss = s.split("\\s+");
                  ContentList.add(ss);
              });
        
        //Generate data for computing.
        P1_Array = ContentList.stream()
                .mapToDouble(s -> Math.log10(Double.parseDouble(s[pv_index[0]])))
                .toArray();
        P2_Array = ContentList.stream()
                .mapToDouble(s -> Math.log10(Double.parseDouble(s[pv_index[1]])))
                .toArray();
        POS_Array = ContentList.stream()
                .mapToInt(s -> Integer.parseInt(s[POS_index]))
                .toArray();
        results = new double[P1_Array.length];
        
        //output title.
        String mm = ContentList.stream()
                .map(s->s[MARKID_COL])
//                .map(Object::toString)
                .collect(Collectors.joining("\t"));
        System.out.println(mm);
        ContentList.clear();
        
      
        setWindowBoundary();
        
        if (Permute_P1 > 0) {
            for (int i = 0; i < Permute_P1; i++) {
                //shuffle array in place
                FisherYatesArrayShuffle.Shuffle(P1_Array);
                outputAllscores();
            }
        }else if (Permute_P2 > 0) {
            for (int i = 0; i < Permute_P2; i++) {
                FisherYatesArrayShuffle.Shuffle(P2_Array);
                outputAllscores();
            }
        }else{
            outputAllscores();
        }
    }
}
