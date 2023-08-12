import java.util.Random;

/**
 * The UncodedBpsk class represents a simulation 
 * of an uncoded Binary Phase-Shift Keying (BPSK) 
 * communication system. 
 * 
 * @author Haiyue Zhang
 */

public class UncodedBpsk {
    public double[] bpsk_ber_sim(int EbNodB_min, int EbNodB_max){
        final int R = 1;
        int n = 1000000;
        
        double[] BER_sim = new double[EbNodB_max - EbNodB_min + 1];
        Random random = new Random();
        double[] msg = new double[n];
        Noise noise = new Noise();
        double[] msg_cap = new double[n];

        int k = 0;

        for (double EbNo_dB = EbNodB_min; EbNo_dB <= EbNodB_max; EbNo_dB++){
            double EbNo = Math.pow(10, (EbNo_dB / 10));
            double sigma = Math.sqrt(1 / (2 * R * EbNo));
            int num_err = 0;
            
            for (int i = 0; i < n; i++){
                msg[i]  = random.nextBoolean() ? 1.0 : 0.0; 
                double m = 1 - 2 * msg[i];
                
                double r = m + sigma * noise.generate_rand_num(1, 0, 1)[0];
                msg_cap[i] = r < 0 ? 1.0 : 0.0;
                num_err += msg[i] == msg_cap[i] ? 0 : 1;
            }
            
            double a = num_err;
            double b = n;
            BER_sim[k] = a/b;
            k++;
        }
        return BER_sim;
    }
}
