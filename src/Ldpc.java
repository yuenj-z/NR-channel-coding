import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

/** 
 * The ldpc class represents a Low-Density Parity-Check (LDPC) 
 * code simulation and encoding/decoding functionality. 
 * LDPC codes are a type of linear error-correcting codes 
 * widely used in modern communication systems to achieve 
 * reliable data transmission over noisy channels.
 * 
 * @author Haiyue Zhang
 */

public class Ldpc {
    Random random = new Random();
    CommFunc func = new CommFunc();
    final int max_itrs = 20; 
    final double s_factor = 0.7; 

    public boolean ldpc_sim(double EbNo_dB, int BG, int Z, int i_LS, double rate){        
        int m; // row
        int n; // col
        
        if (BG == 1){
            m = 46;
            n = 68;    
        } 
        else {
            m = 42;
            n = 52;     
        }

        int k = n - m; 
        int nbRM = (int)(k / rate) + 2; 
        int mbRM = Math.max(4, nbRM - k); 

        /* ENCODER */
        double[] cword = new double[n * Z]; 
        double[] RMcword = new double[nbRM * Z]; 
        double[] msg = new double[k * Z]; 
       
        for (int i = 0; i < k * Z; i++) {
            msg[i]  = random.nextBoolean() ? 1.0 : 0.0; 
        } 
        
        cword = encoder(Z, i_LS, msg, BG); 
        RMcword = Arrays.copyOfRange(cword, 0, nbRM * Z - 1); 
        
        /* MODULATION, CHANNEL */
        Noise noise = new Noise();
        double[] r = new double[nbRM * Z]; 
        double EbNo = Math.pow(10, (EbNo_dB/10)); 
        double sigma = Math.sqrt(1/(2 * k / (nbRM-2) * EbNo)); 

        for (int i = 0; i < RMcword.length; i++){ 
            double mod = 1 - 2 * RMcword[i]; 
            r[i] = mod + sigma * noise.generate_rand_num(1, 0, 1)[0];      
        }   
        
        for (int i = 0; i < 2*Z; i++){
            r[i] = 0; 
        }
        
        /* DECODER */
        return decoder(r , BG, Z, i_LS, rate, mbRM, k, msg);
    }

    public double[] encoder(int Z, int i_LS, double[] msg, int BG){
        int m; 
        int n; 
        double[] temp = new double[Z]; 
        HashMap<Integer, int[]> H; 

        if (BG == 1){
            m = 46;
            n = 68;
            H = LdpcTables.H_bg1;
        } 
        else {
            m = 42;
            n = 52;
            H = LdpcTables.H_bg2; 
        }

        double[] cword = new double[n * Z];
        
        System.arraycopy(msg, 0, cword, 0, msg.length); 

        /* GENERATE PARITY BITS */
        for (int i = 0; i < 4; i++) { 
            for (int j = 0; j <= n - m - 1; j++){
                temp = hct(H, Z, i, j, cword, temp, i_LS); 
            }
        } 

        int i = 1;
        int j = n - m;
        int p0_sh;

        int[] default_ = {-1};
        
        if (H.getOrDefault(i * 100 + j, default_)[0] == -1){ 
            p0_sh = H.get((i + 1) * 100 + j)[i_LS];
        } 
        else {
            p0_sh = H.get(i * 100 + j)[i_LS];
        }

        System.arraycopy(mul_sh(temp, (Z - p0_sh % Z) % Z ), 0, cword, (n - m) * Z, Z); 

        for (i = 0; i < 3; i++) { 
            for (int z = 0; z < Z; z++){
                temp[z] = 0; 
            }
            for (int y = 0; y <= n - m + i; y++){
                temp = hct(H, Z, i, y, cword, temp, i_LS); 
            }
            System.arraycopy(temp, 0, cword, (n - m + i + 1) * Z, Z);
        } 

        for (i = 4; i < m; i++) { 
            for (int z = 0; z < Z; z++){
                temp[z] = 0;
            }
            for (int y = 0; y < n - m + 4; y++){
                temp = hct(H, Z, i, y, cword, temp, i_LS);
            }
            System.arraycopy(temp, 0, cword, (n - m + i) * Z, Z);
        } 
        return cword;
    }

    public double[] mul_sh(double[] m, int i){
        if (i < 0 && i != -1){ 
            throw new UnsupportedOperationException();
        }

        double[] y = new double[m.length];

        if (i == -1){
            for (int x = 0; x < m.length; x++){
                y[x] = 0; 
            }
        } 
        else {
            for (int x = 0; x < m.length; x++){
                y[x] = m[i];
                if (i == m.length - 1){
                    i = 0;
                } 
                else{
                    i++;
                } 
            }
        }
        return y;
    }

    public double[] hct(HashMap<Integer, int[]> H, int Z, int i, int j, double[] cword, double[] temp, int i_LS){
        int key = i * 100 + j;
        int v_i_j = -1; 
        
        if (H.get(key) != null){
            v_i_j = H.get(key)[i_LS];
        }
        
        int shift = v_i_j % Z;

        double[] t1 = mul_sh(Arrays.copyOfRange(cword, j * Z, (j + 1) * Z), shift);

        temp = func.add(temp, t1);
   
        for (int x = 0; x < temp.length; x++) {
            temp[x] = temp[x] % 2.0;
        } 
        return temp;
    }

    public boolean decoder(double[] r, int BG, int Z, int i_LS, double rate, int mbRM, int k, double[] msg){
        double[] L = new double[r.length]; 
        double[][] R = new double[totCol(BG)][Z]; 
        double[][] Treg = new double[maxCol(BG)][Z]; 
        HashMap<Integer, int[]> H;
        int[][] bg;

        System.arraycopy(r, 0, L, 0, r.length);

        if (BG == 1){
            bg = LdpcTables.col_bg1;
            H = LdpcTables.H_bg1;
        }
        else{
            bg = LdpcTables.col_bg2;
            H = LdpcTables.H_bg2;
        }

        double[] temp = new double[Z];
        
        for (int itr = 1; itr <= max_itrs; itr++){
            int tcol = 0;

            for (int i = 0; i < mbRM; i++){
                int rcol = 0;

                for (int j : bg[i]) { 
                    temp = Arrays.copyOfRange(L, j * Z, (j + 1) * Z); 
                    temp = func.sub(temp, R[tcol]);
                    System.arraycopy(temp, 0, L, j * Z, Z);

                    int key = i * 100 + j;
                    int v_i_j = -1; 
        
                    if (H.get(key) != null){
                        v_i_j = H.get(key)[i_LS];
                    }

                    Treg[rcol] = mul_sh(temp, v_i_j % Z); 
                    
                    tcol++; 
                    rcol++; 
                }
                
                /* MINSUM */
                double[] sign = new double[rcol];
                double[] value = new double[rcol]; 

                for (int z = 0; z < Z; z++){
                    double prod_s = 1.0;
                    for (int d = 0; d < rcol; d++) {
                        
                        sign[d] = func.get_sign(Treg[d][z]);
                        value[d] = Treg[d][z];                        
                        prod_s *= sign[d];
                    }
                    
                    int m1_index = min(value, -1);
                    double m1 = Math.abs(value[m1_index]);
                    int m2_index = min(value, m1_index);
                    double m2 = Math.abs(value[m2_index]);
                                      
                    for (int d = 0; d < rcol; d++){
                        if (d == m1_index){
                            Treg[d][z] = m2 * prod_s * sign[d] * s_factor;
                        }
                        else{
                            Treg[d][z] = m1 * prod_s * sign[d] * s_factor;
                        }                                               
                    }                    
                } 

                tcol -= rcol;
                rcol = 0;

                for (int j : bg[i]) {
                    int key = i * 100 + j;
                    int v_i_j = -1; 
        
                    if (H.get(key) != null){
                        v_i_j = H.get(key)[i_LS];
                    }
                    for (int z = 0; z < Z; z++){
                        temp[z] = Treg[rcol][z];
                    }
                    if (v_i_j% Z == 0){
                        R[tcol] = mul_sh(temp, 0);
                    }
                    else{
                        R[tcol] = mul_sh(temp, Z - v_i_j%Z); 
                    }

                    temp = Arrays.copyOfRange(L, j*Z, (j+1)*Z); 
                    temp = func.add(temp, R[tcol]);
                    System.arraycopy(temp, 0, L, j*Z, Z);  
                    tcol++;
                    rcol++;                  
                }
            }

            double[] msg_cap = new double[k * Z];
           
            for (int c = 0; c < msg_cap.length; c++){
                msg_cap[c] = L[c] < 0 ? 1.0 : 0.0;               
            }

            int sum = 0;
            for (int c = 0; c < msg_cap.length; c++){
                sum += Math.abs(msg_cap[c] - msg[c]);
            }

            if (sum == 0){
                return true;
            }         
        }
        return false;
    }

    public void test_ldpc_encode(int i_LS, int BG){
        Ldpc ldpc = new Ldpc();
        LdpcTables ldpc_tables = new LdpcTables();
        
        int m; 
        int n; 
        
        if (BG == 1) { 
            n = 68; 
            m = 46; 
        } 
        else { 
            n = 52; 
            m = 42; 
        }
        
        int k = n-m;
        
        for (int Z : ldpc_tables.zset[i_LS]) {
            double[] cword = new double[n * Z];
            double[] msg = new double[k * Z];
            
            for (int i = 0; i < k * Z; i++) {
                msg[i]  = random.nextBoolean() ? 1.0 : 0.0;                
            }
            
            cword = ldpc.encoder(Z, i_LS, msg, BG);
        }
    }

    public boolean checkcword(int Z, double[] cword, int BG, int i_LS){
        int m;
        int n;
        HashMap<Integer, int[]> H;

        if (BG == 1){
            m = 46;
            n = 68;
            H = LdpcTables.H_bg1;
        } else {
            m = 42;
            n = 52;
            H = LdpcTables.H_bg2; 
        }

        double[] syn = new double[m * Z];
        double[] temp = new double[Z];

        for (int i = 0; i < m; i++){
            for (int z = 0; z < Z; z++){
                temp[z] = 0;
            }
            for (int j = 0; j < n; j++){
                temp = hct(H, Z, i, j, cword, temp, i_LS);   
            }
            System.arraycopy(temp, 0, syn, i * Z , Z);
        }

        double synsum = 0.0;
        for (int i = 0; i < syn.length; i++){
            synsum = synsum + syn[i];
        }
        return synsum == 0;
    }

    public int totCol(int BG){
        int[][] bg;

        if (BG == 1){
            bg = LdpcTables.col_bg1;
        }
        else{
            bg = LdpcTables.col_bg2;
        }
           
        int count = 0;

        for (int i = 0; i < bg.length; i++){
            count += bg[i].length;
        }
        return count;
    }

    public int maxCol(int BG){
        int[][] bg;

        if (BG == 1){
            bg = LdpcTables.col_bg1;
        }
        else{
            bg = LdpcTables.col_bg2;
        }
           
        int max = 0;
        
        for (int i = 0; i < bg.length; i++){
            if (bg[i].length > max){
                max = bg[i].length;
            }
        }
        return max;
    }
    
    public int min(double[] arr, int m1_index){
        int index = 0;
        double min = 1000;
        
        for (int i = 0; i < arr.length; i++){
            if (i != m1_index && Math.abs(arr[i]) < min){
                min = Math.abs(arr[i]);
                index = i;
            }
        }
        return index;
    }

    public static void main(String[] args) {
        Ldpc ldpc = new Ldpc();
        LdpcTables ldpc_tables = new LdpcTables();

        double[] ebno = {1, 1.5, 1.75, 2};

        /***
         * CHANGE YOUR VARIABLES HERE ------------------------->
         */

        int nsim = 2000; // number of simulation to run  
        int BG = 2; // choose 1 or 2
        int Z = 52; 
        int i_LS = 6; 
        double rate = 0.5;
        
        /* Pick your i_LS based on your Z
         * Z =  { 2,  4,  8,  16, 32, 64, 128, 256 }, //i_LS = 0
                { 3,  6, 12,  24, 48, 96, 192, 384 }, //i_LS = 1
                { 5, 10, 20,  40, 80, 160, 320},      //i_LS = 2
                { 7, 14, 28,  56, 112, 224},          //i_LS = 3
                { 9, 18, 36,  72, 144, 288},          //i_LS = 4
                {11, 22, 44,  88, 176, 352},          //i_LS = 5
                {13 , 26, 52, 104, 208},              //i_LS = 6
                {15, 30, 60, 120, 240}                //i_LS = 7
        */

        /**
         * <---------------------------------------------------
         */

        for (double d : ebno) {
            double num_not_ok = 0;
            for (int i = 0; i < nsim; i++){
                if (!ldpc.ldpc_sim(d, BG, Z , i_LS, rate)){
                    num_not_ok++;
                }
            }
            System.out.println("bler for " + d + ": " + num_not_ok/2000.0);
        }
    }
}


