public class BerTh {
    public double[] bpsk_ber_th(int EbNodB_min, int EbNodB_max){
        final int R = 1;
        double[] BER_th = new double[EbNodB_max - EbNodB_min + 1];

        int i = 0;

        for (double EbNo_dB = EbNodB_min; EbNo_dB <= EbNodB_max; EbNo_dB++){
            double EbNo = Math.pow(10, (EbNo_dB / 10));
            double sigma = Math.sqrt(1/(2 * R * EbNo));

            BER_th[i] = 0.5 * (ErrorFunction.erfc(Math.sqrt(EbNo)));
            
            i++;
        }
        return BER_th;
    }
}
