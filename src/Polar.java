import java.lang.Math;
import java.util.Random;

public class Polar {
    Random random = new Random();
    CommFunc func = new CommFunc();

    public boolean polar_sim(int N, double EbNo_dB, int K){
        /* ENCODER */
        double[] msg = new double[K];
        double[] u = new double[N];
        for (int i = 0; i < K; i++) {
            msg[i]  = random.nextBoolean() ? 1.0 : 0.0; 
        }
        u = encoder(N, K, msg);

        /* MODULATION, CHANNEL */
        Noise noise = new Noise();
        double[] r = new double[N]; 
        double EbNo = Math.pow(10, (EbNo_dB/10)); 
        double sigma = Math.sqrt(1/(2 * ((double)K/(double)N) * EbNo)); 

        for (int i = 0; i < N; i++){ 
            double mod = 1 - 2 * u[i]; 
            r[i] = mod + sigma * noise.generate_rand_num(1, 0, 1)[0];      
        }   

        /* DECODER */
        return decoder(N, r, K, msg);
    }

    public double[] encoder(int N, int K, double[] msg){
        double n = Math.log(N) / Math.log(2);
        if(!(n == Math.floor(n))){
            throw new UnsupportedOperationException("N has to be 2^n");
        }
        
        double[] u = new double[N];
        int i = 0;
        int[] q = Q(N);

        for(int j = N - K; j < N; j++){
            u[q[j]] = msg[i];
            i++;
        }     

        int m = 1;
        
        for (int d = (int)n - 1; d >= 0; d--){
            for (int j = 0; j < N; j = j + 2 * m){
            
                double[] a = new double[m];
                double[] b = new double[m];
                double[] temp = new double[m];
             
                System.arraycopy(u, j, a, 0, m);
                System.arraycopy(u, j + m, b, 0, m);
                temp = func.add(a, b);

                for (int x = 0; x < temp.length; x++) {
                    temp[x] = temp[x] % 2.0;
                }

                System.arraycopy(temp, 0, u, j, m);
            }
            m *= 2;
        }  
        return u;
    }

    public int[] Q(int N){
        PolarTables polar_tables = new PolarTables();
        int[] Q = new int[N];
        int j = 0;

        for(int i = 0; i < polar_tables.Qi.length; i++){
            if (polar_tables.Qi[i] < N){
                Q[j] = polar_tables.Qi[i];
                j++;
            }
        }
        return Q;
    }

    public boolean decoder(int N, double[] r, int K, double[] msg){
        int n = (int)(Math.log(N) / Math.log(2));
        if (Math.pow(2, n) != N){
            throw new UnsupportedOperationException("N has to be 2^n");
        }
        double[][] L = new double[n + 1][N];
        double[][] ucap = new double[n + 1][N];
        int[] ns = new int[2 * N - 1]; 

        for (int i = 0; i < N; i++){
            L[0][i] = r[i];
        }

        int node = 0;
        int depth = 0;
        boolean done = false;

        int[] Q = Q(N);
        int[] F = new int[N - K];
        System.arraycopy(Q, 0, F, 0, N - K);

        while (!done){

            /* LEAF NODE */
            if (depth == n){  
                if (func.contains(F, node)){
                    ucap[n][node] = 0;
                }
                else{
                    ucap[n][node] = L[n][node] >= 0 ? 0 : 1;
                }
                if (node == N - 1){
                    done = true;
                }
                else{
                    node = node / 2;
                    depth--;
                }
            }

            /* NOT LEAF NODE */
            else{ 
                int npos = (int)Math.pow(2, depth) + node - 1;
                int temp = (int)Math.pow(2, n - depth);

                /* STEP L */
                if (ns[npos] == 0){ 
                    double[] ln = new double[temp];
                    System.arraycopy(L[depth], temp * node, ln, 0, temp);
                    double[] a = splitL(ln);
                    double[] b = splitR(ln);

                    node *= 2;
                    depth++;
                    temp /= 2;

                    System.arraycopy(f(ln), 0, L[depth], temp * node, temp);

                    ns[npos] = 1;
                }

                /* STEP R */
                else if(ns[npos] == 1){ 
                    double[] ln = new double[temp];
                    System.arraycopy(L[depth], temp * node, ln, 0, temp);
                    double[] a = splitL(ln);
                    double[] b = splitR(ln);

                    int lnode = 2 * node;
                    int ldepth = depth + 1;
                    int ltemp = temp / 2;

                    double[] ucapn = new double[ltemp];
                    System.arraycopy(ucap[ldepth], ltemp * lnode, ucapn, 0, ltemp);

                    node = node * 2 + 1;
                    depth++;
                    temp = temp / 2;

                    System.arraycopy(g(ln, ucapn), 0, L[depth], temp * node, temp);
    
                    ns[npos] = 2;
                }

                /* STEP U */
                else if(ns[npos] == 2){ 
                    int lnode = 2 * node;
                    int rnode = 2 * node + 1;
                    int cdepth = depth + 1;
                    int ctemp = temp / 2; 

                    double[] ucapl = new double[ctemp];
                    double[] ucapr = new double[ctemp];

                    System.arraycopy(ucap[cdepth], ctemp * lnode, ucapl, 0, ctemp);
                    System.arraycopy(ucap[cdepth], ctemp * rnode, ucapr, 0, ctemp);
                    double[] ucapn = func.add(ucapl, ucapr);
                    for (int x = 0; x < ucapn.length; x++) {
                        ucapn[x] = ucapn[x] % 2.0;
                    }

                    System.arraycopy(ucapn, 0, ucap[depth], temp * node, ctemp);
                    
                    
                    System.arraycopy(ucapr, 0, ucap[depth], temp * node + ctemp, ctemp);

                    node = node / 2; 
                    depth--;
                }
            }
        }

        double[] msg_cap = new double[K];

        int i = 0;
        for(int j = N - K; j < N; j++){
            msg_cap[i] = ucap[n][Q[j]];
            i++;   
        }

        int sum = 0;

        for (int c = 0; c < msg_cap.length; c++){
            sum += Math.abs(msg_cap[c] - msg[c]);
        }
        if (sum == 0){
            return true;
        }
        return false;
    }

    public double[] splitL(double[] u){
        if (u.length % 2 != 0){
            throw new UnsupportedOperationException("u has to be an even number");
        }
        int length = u.length / 2;
        double[] a = new double[length];
        
        System.arraycopy(u, 0, a, 0, length);
        return a;
    }

    public double[] splitR(double[] u){
        if (u.length % 2 != 0){
            throw new UnsupportedOperationException("u has to be an even number");
        }
        int length = u.length / 2;
        double[] b = new double[length];
        
        System.arraycopy(u, length, b, 0, length);
        return b;
    }

    public double[] f(double[] u){
        double[] a = splitL(u);
        double[] b = splitR(u);
        double[] f = new double[a.length];

        for (int i = 0; i < a.length; i++){
            double sign = func.get_sign(a[i]) * func.get_sign(b[i]);
            if (Math.abs(a[i]) < Math.abs(b[i])){
                f[i] = Math.abs(a[i]) * sign;
            }
            else {
                f[i] = Math.abs(b[i]) * sign;
            }
        }
        return f;
    }

    public double[] g(double[] u, double[] c){
        double[] a = splitL(u);
        double[] b = splitR(u);
        
        double[] g = new double[a.length];

        if (a.length != c.length){
            throw new UnsupportedOperationException("a has to be same length as c");
        }

        for (int i = 0; i < a.length; i++){
            g[i] = b[i] + (1 - 2 * c[i]) * a[i]; 
        }
        return g;
    }
     

    public static void main(String[] args) {
        Polar polar = new Polar();

        double[] ebnoList = {1, 1.5, 2.0, 2.5, 3.0};

        /***
         * CHANGE YOUR VARIABLES HERE ------------------------->
         */

        int nsim = 2000; // number of simulation to run  
        int N = 1024;
        int K = 512;

        /**
         * <----------------------------------------------------
         */

        for (double ebno : ebnoList) {
            double num_not_ok = 0;
            for (int i = 0; i < nsim; i++){
                if (!polar.polar_sim(N, ebno, K)){
                    num_not_ok++;
                }
            }
            System.out.println("bler for " + ebno + ": " + num_not_ok/(double)nsim);
        }
        
    }
}
